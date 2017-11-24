#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use SAM_reader;
use SAM_entry;
use DelimParser;
use Carp;
use Data::Dumper;
use Set::IntervalTree;
use SeqUtil;
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);



my $MIN_ALIGN_PER_ID = 96;

my $MIN_SMALL_ANCHOR = 12;
my $MIN_LARGE_ANCHOR = 25;
my $MAX_END_CLIP = 10;
my $MIN_SEQ_ENTROPY = 1.5;

my $DEBUG = 0;

my $usage = <<__EOUSAGE__;

###############################################################
#
# Required:
#
#  --gtf_file <string>     genePairContig.gtf
#  --bam <string>          read_alignments.bam
#
# Optional:
#
#  --MIN_ALIGN_PER_ID <int>   default: $MIN_ALIGN_PER_ID
#
#  --MIN_SMALL_ANCHOR <int>   default: $MIN_SMALL_ANCHOR
#  --MIN_LARGE_ANCHOR <int>   default: $MIN_LARGE_ANCHOR
#  --MAX_END_CLIP <int>       default: $MAX_END_CLIP
#  --MIN_SEQ_ENTROPY <float>    default: $MIN_SEQ_ENTROPY
#
#  -d                         debug mode
#
##############################################################


__EOUSAGE__

    ;


my $gtf_file;
my $bam_file;
my $help_flag;

&GetOptions('help|h' => \$help_flag,
            'gtf_file=s' => \$gtf_file,
            'bam=s' => \$bam_file,
            
            'MIN_ALIGN_PER_ID=i' => \$MIN_ALIGN_PER_ID,
            
            'MAX_END_CLIP=i' => \$MAX_END_CLIP,
            'MIN_SMALL_ANCHOR=i' => \$MIN_SMALL_ANCHOR,
            'MIN_LARGE_ANCHOR=i' => \$MIN_LARGE_ANCHOR,
            'MIN_SEQ_ENTROPY=f' => \$MIN_SEQ_ENTROPY,

            'd' => \$DEBUG,
  );


if ($help_flag) {
    die $usage;
}

unless ($gtf_file && $bam_file) {
    die $usage;
}


main: {
    
    my %scaffold_itree_to_exon_structs;
    print STDERR "-parsing GTF file: $gtf_file\n";
    my %scaffold_to_gene_structs = &parse_gtf_file($gtf_file, \%scaffold_itree_to_exon_structs);
        
    print STDERR Dumper(\%scaffold_to_gene_structs) if $DEBUG;
    
    my $prev_scaff = "";

    my ($left_gene, $right_gene, $gene_bound_left, $gene_bound_right);
    
    my %fusion_split_reads;
    
    my %fusion_junctions;
    my %fusion_large_anchors;

    my %sam_index_capture;
    
    my $counter = 0;
    ## find the reads that matter:
    my $sam_reader = new SAM_reader($bam_file);
    while (my $sam_entry = $sam_reader->get_next()) {

        if ($DEBUG) {
            print $sam_entry->get_original_line() . "\n";
        }
        
        $counter++;
        if ($counter % 10000 == 0) { 
            print STDERR "\r[$counter]   ";
        }

        if ($sam_entry->is_duplicate()) {
            
            if ($DEBUG) { print STDERR "-skipping duplicate mapping\n"; }
            next;
        }
                
        my $qual_val = $sam_entry->get_mapping_quality();
        unless ($qual_val > 0) { 
            if ($DEBUG) { print STDERR "-skipping due to mapping qual <= 0\n"; }
            next; 
        }
        

        ## examine number of mismatches in read alignment
        my $mismatch_count = 0;
        my $line = $sam_entry->get_original_line();

        if (0) {  # turning off - allow for multimapping reads
            ## ensure the match is unique
            $line =~ /NH:i:(\d+)/ or die "Error, cannot extract hit count (NH:i:) from sam entry: $line";
            my $num_hits = $1;
            if ($num_hits != 1) {
                if ($DEBUG) { print STDERR "-skipping, num hits ($num_hits) indicates not unique\n"; }
                next;
            }
        }
        
        if ($line =~ /NM:i:(\d+)/i) {
            $mismatch_count = $1;
        }
        my $alignment_length = $sam_entry->get_alignment_length();
        unless ($alignment_length) {
            if ($DEBUG) { print STDERR "-skipping, no alignment length\n"; }
            next;
        }
        my $per_id = ($alignment_length - $mismatch_count) / $alignment_length * 100;
        if ($per_id < $MIN_ALIGN_PER_ID) {
            if ($DEBUG) { print STDERR "-skipping, per_id $per_id < min required: $MIN_ALIGN_PER_ID\n";}
            next;
        }
        
        ## check end clipping of alignment
        my $cigar = $sam_entry->get_cigar_alignment();
        if (&alignment_has_excessive_soft_clipping($cigar)) {
            if ($DEBUG) { print STDERR "-skipping, excessive soft clipping ($cigar)\n";  }
            next;
        }
        
        my $read_name = $sam_entry->reconstruct_full_read_name();
        
        my $core_read_name = $sam_entry->get_core_read_name();
        
        my $read_sequence = $sam_entry->get_sequence();
        
        my $scaffold = $sam_entry->get_scaffold_name();
        if ($scaffold ne $prev_scaff) {
            
            if (! exists $scaffold_to_gene_structs{$scaffold}) { 
                if ($DEBUG) { print STDERR "-skipping, not a known fusion contig target\n"; }
                next;  # STAR aligns to the whole genome + fusion contigs.
            }
            my @gene_structs = @{$scaffold_to_gene_structs{$scaffold}};
            
            @gene_structs = sort {$a->{lend} <=> $b->{lend}} @gene_structs;

            if (scalar @gene_structs > 2) {
                @gene_structs = &partition_gene_structs(@gene_structs);
            }
            
            
            
            $left_gene = $gene_structs[0];
            $right_gene = $gene_structs[1];
            
            $gene_bound_left = $left_gene->{rend};
            $gene_bound_right = $right_gene->{lend};
            
            
            if ($gene_bound_left > $gene_bound_right) {
                die "Error, bounds out of order: $gene_bound_left ! <  $gene_bound_right";
            }
            
            
            $prev_scaff = $scaffold;
        }
        
        my ($span_lend, $span_rend) = $sam_entry->get_genome_span();
        if ($span_lend < $gene_bound_left && $span_rend > $gene_bound_right) {
            
            my ($genome_coords_aref, $read_coords_aref) = $sam_entry->get_alignment_coords();
            
            my $itree = $scaffold_itree_to_exon_structs{$scaffold};
            my @exon_hits = &hits_exon_bound($genome_coords_aref, $read_coords_aref, $itree);
            
            print Dumper(\@exon_hits) if $DEBUG;

            unless (scalar @exon_hits > 1) { next; } 
            
            my %genes_matched;
            foreach my $exon_hit (@exon_hits) {

                my $gene_id = $exon_hit->{exon_struct}->{gene_id};
                $gene_id =~ s/\^.*//;
                $genes_matched{ $gene_id }++;
            }
            my $num_genes_matched = scalar(keys %genes_matched);
            print STDERR "Genes matched: $num_genes_matched\n" if $DEBUG;
            unless ($num_genes_matched > 1) { next; }

            print STDERR "Genes matched: " . join(" ", keys %genes_matched) . "\n" if $DEBUG;
            
            my $junction_coord_token = &get_junction_coord_token(@exon_hits);
            if ($junction_coord_token) {
                
                my ($gene_A, $brkpt_A, $gene_B, $brkpt_B, $splice_token) = split(/\t/, $junction_coord_token);
                my ($left_read_segments_aref, $right_read_segments_aref) = &divide_junction_read_at_breakpoint($genome_coords_aref, $read_coords_aref, $brkpt_A);
                
                # coordinate sets based on read as reference coordinates
                my @left_read_coordsets = sort {$a->[0]<=>$b->[0]} @$left_read_segments_aref;
                my @right_read_coordsets = sort {$a->[0]<=>$b->[0]} @$right_read_segments_aref;

                my $left_immediate_breakpoint_segment = pop @left_read_coordsets;
                my $right_immediate_breakpoint_segment = shift @right_read_coordsets;

                my $left_brkpt_length = abs($left_immediate_breakpoint_segment->[1] - $left_immediate_breakpoint_segment->[0]) + 1;
                my $right_brkpt_length = abs($right_immediate_breakpoint_segment->[1] - $right_immediate_breakpoint_segment->[0]) + 1;
                
                print STDERR "LEFT_read_seg: $left_brkpt_length, RIGHT_read_seg: $right_brkpt_length\n" if $DEBUG;

                # anchor length check
                unless ($left_brkpt_length >= $MIN_SMALL_ANCHOR
                    &&
                    $right_brkpt_length >= $MIN_SMALL_ANCHOR) {
                    next;
                }

                ## anchor sequence entropy check
                unless (&has_min_anchor_seq_entropy(\$read_sequence, $left_immediate_breakpoint_segment, $MIN_SEQ_ENTROPY)
                        &&
                        &has_min_anchor_seq_entropy(\$read_sequence, $right_immediate_breakpoint_segment, $MIN_SEQ_ENTROPY) ) {
                    
                    next;
                }
                
                # calling it a fusion read.
                $fusion_split_reads{$core_read_name} = 1;
                $sam_index_capture{$counter} = 1;
                push (@{$fusion_junctions{$junction_coord_token}}, $read_name);
                
                
                if ($left_brkpt_length >= $MIN_LARGE_ANCHOR) {
                    $fusion_large_anchors{$junction_coord_token}->{LEFT_LARGE_ANCHOR}++;
                }
                if ($right_brkpt_length >= $MIN_LARGE_ANCHOR) {
                    $fusion_large_anchors{$junction_coord_token}->{RIGHT_LARGE_ANCHOR}++;
                }
                
                
                
            }
        }
        

    }
    

    {
        # capture the alignments involving identified junction / split reads:

        my $sam_reader = new SAM_reader($bam_file);
        my $counter = 0;
        while (my $sam_entry = $sam_reader->get_next()) {

            $counter++;
            
            my $core_read_name = $sam_entry->get_core_read_name();
            # want both reads in the pair
            if ($fusion_split_reads{$core_read_name} && $sam_index_capture{$counter}) {
                print $sam_entry->get_original_line() . "\n";
            }
        }
    
    }
    


    {
        
        print STDERR "-writing fusion junction support info.\n";
        
        ## output the junction info summary:
        
        my $junction_coords_info_file = "$bam_file.fusion_junction_info";
        open (my $ofh, ">$junction_coords_info_file") or die "Error, cannot write to file: $junction_coords_info_file";
        
        my @column_headers = qw(LeftGene LeftLocalBreakpoint LeftBreakpoint
                                RightGene RightLocalBreakpoint RightBreakpoint
                                SpliceType JunctionReadCount LargeAnchorSupport JunctionReads);
        
        my $tab_writer = new DelimParser::Writer($ofh, "\t", \@column_headers);
        
        my @fusion_j = reverse sort { $#{$fusion_junctions{$a}} <=> $#{$fusion_junctions{$b}} } keys %fusion_junctions;
        foreach my $fusion (@fusion_j) {
            my @reads = @{$fusion_junctions{$fusion}};
            my $num_fusion_reads = scalar @reads;
            
            my $LARGE_ANCHOR_BOTH_ENDS = ($fusion_large_anchors{$fusion}->{LEFT_LARGE_ANCHOR}
                                          &&
                                          $fusion_large_anchors{$fusion}->{RIGHT_LARGE_ANCHOR})
                ? "YES"
                : "NO";
            
            my ($LeftGene, $LeftLocalBreakpoint, $LeftBreakpoint,
                $RightGene, $RightLocalBreakpoint, $RightBreakpoint,
                $Splice_type) = split(/\t/, $fusion);
            
            $tab_writer->write_row( { LeftGene => $LeftGene,
                                      LeftLocalBreakpoint => $LeftLocalBreakpoint,
                                      LeftBreakpoint => $LeftBreakpoint,
                                      RightGene => $RightGene,
                                      RightLocalBreakpoint => $RightLocalBreakpoint,
                                      RightBreakpoint => $RightBreakpoint,
                                      SpliceType => $Splice_type,
                                      JunctionReadCount => $num_fusion_reads,
                                      LargeAnchorSupport => $LARGE_ANCHOR_BOTH_ENDS,
                                      JunctionReads => join(",", @reads) } );
            
            
            
            
        }
        
        close $ofh;
    }
    


    exit(0);
}



####
sub get_junction_coord_token {
    my @exon_hits = @_;
    
    @exon_hits = sort {$a->{align_genome_lend}<=>$b->{align_genome_lend}} @exon_hits;
    
    print STDERR "Ordered exon hits: " . Dumper(\@exon_hits) if $DEBUG;

    
    my @fusion_candidates;
    
    # find the pair that fuses the gene
    for (my $i = 0; $i < $#exon_hits; $i++) {
        
        my $exon_hit_i = $exon_hits[$i];
        my $exon_struct_i = $exon_hit_i->{exon_struct};

        my $gene_id_i = $exon_struct_i->{gene_id};
        $gene_id_i =~ s/\^.*//;
        
        for (my $j = $i+1; $j <= $#exon_hits; $j++) {
            
            my $exon_hit_j = $exon_hits[$j];
            my $exon_struct_j = $exon_hit_j->{exon_struct};

            my $gene_id_j = $exon_struct_j->{gene_id};
            $gene_id_j =~ s/\^.*//;
            
            unless ($exon_struct_i->{rend} < $exon_struct_j->{lend}) { next; }
            
            ## see if different genes and properly split alignment
            if ($gene_id_i ne $gene_id_j
                &&
                ( abs($exon_hit_i->{read_end3} - $exon_hit_j->{read_end5}) == 1 # forward read alignment
                  ||
                  abs($exon_hit_i->{read_end5} - $exon_hit_j->{read_end3}) == 1 # reverse read alignment
                ) 
                ) {
                
                ## determine deltas between exon ends and read alignment ends
                my $delta_i = $exon_hit_i->{align_genome_rend} - $exon_hit_i->{exon_struct}->{rend};
                my $delta_j = $exon_hit_j->{align_genome_lend} - $exon_hit_j->{exon_struct}->{lend};
                my $score = sqrt( $delta_i**2 + $delta_j ** 2);
                
                push (@fusion_candidates, {
                    
                    exon_hit_i => $exon_hit_i, 
                    exon_hit_j => $exon_hit_j,

                    delta_i => $delta_i,
                    delta_j => $delta_j,

                    score => $score,
                    
                    });
            }
        }
    }
    

    unless (@fusion_candidates) {
        print STDERR "** no fusion candidate.\n" if $DEBUG;
        return(undef);
    }


    @fusion_candidates = sort {$a->{score} <=> $b->{score}} @fusion_candidates;


    my $best_fusion_candidate = shift @fusion_candidates;
    
    print Dumper($best_fusion_candidate) if $DEBUG;
    
    my $splice_token = ($best_fusion_candidate->{score} == 0) ? "ONLY_REF_SPLICE" : "INCL_NON_REF_SPLICE";


    ## determine genome coordinate based on fusion-contig coordinate for breakpoints
    
    my $left_breakpoint_coordinate = ($best_fusion_candidate->{exon_hit_i}->{exon_struct}->{orig_orient} eq '+')
        ? ($best_fusion_candidate->{exon_hit_i}->{exon_struct}->{orig_end3} + $best_fusion_candidate->{delta_i})
        : ($best_fusion_candidate->{exon_hit_i}->{exon_struct}->{orig_end3} - $best_fusion_candidate->{delta_i});

    $left_breakpoint_coordinate = join(":", 
                                       $best_fusion_candidate->{exon_hit_i}->{exon_struct}->{orig_chr}, 
                                       $left_breakpoint_coordinate,
                                       $best_fusion_candidate->{exon_hit_i}->{exon_struct}->{orig_orient}); 
    
    my $right_breakpoint_coordinate = ($best_fusion_candidate->{exon_hit_j}->{exon_struct}->{orig_orient} eq '+')
        ? ($best_fusion_candidate->{exon_hit_j}->{exon_struct}->{orig_end5} + $best_fusion_candidate->{delta_j})
        : ($best_fusion_candidate->{exon_hit_j}->{exon_struct}->{orig_end5} - $best_fusion_candidate->{delta_j});
    
    $right_breakpoint_coordinate = join(":", 
                                        $best_fusion_candidate->{exon_hit_j}->{exon_struct}->{orig_chr}, 
                                        $right_breakpoint_coordinate,
                                        $best_fusion_candidate->{exon_hit_j}->{exon_struct}->{orig_orient}); 
    
    

    my $fusion_token = join("\t", 
                            $best_fusion_candidate->{exon_hit_i}->{exon_struct}->{gene_id},
                            #$best_fusion_candidate->{exon_hit_i}->{rend},
                            $best_fusion_candidate->{exon_hit_i}->{align_genome_rend},
                            #$best_fusion_candidate->{exon_hit_i}->{exon_struct}->{orig_end3},
                            $left_breakpoint_coordinate,


                            $best_fusion_candidate->{exon_hit_j}->{exon_struct}->{gene_id},
                            #$best_fusion_candidate->{exon_hit_j}->{lend},
                            $best_fusion_candidate->{exon_hit_j}->{align_genome_lend},
                            #$best_fusion_candidate->{exon_hit_j}->{exon_struct}->{orig_end5},
                            $right_breakpoint_coordinate,
                            
                            $splice_token,

                            #### TODO: determine reference genome coordinates for the fusion coordinate.
        );
    
    return($fusion_token);
    
    
}




####
sub hits_exon_bound {
    my ($genome_coords_aref, $read_coords_aref, $itree) = @_;

    my @read_coordsets = @$read_coords_aref;
    
    my @hits;

    foreach my $coordset (@$genome_coords_aref) {
        my ($lend, $rend) = @$coordset;

        if ($lend == $rend) { next; } # no single base searches

        my $read_coords_aref = shift @read_coordsets;
        my ($read_end5, $read_end3) = @$read_coords_aref;
        
        my $overlaps_aref = $itree->fetch($lend, $rend);
        
        foreach my $overlap_struct (@$overlaps_aref) {
            push (@hits, { exon_struct => $overlap_struct,
                           
                           align_genome_lend => $lend,
                           align_genome_rend => $rend,
                           
                           read_end5 => $read_end5,
                           read_end3 => $read_end3,
                  } );
        }
    }
    
    return(@hits);
    
}


####
sub parse_gtf_file {
    my ($gtf_file, $scaffold_itree_to_exon_structs_href) = @_;

    my %scaff_to_gene_to_coords;
    
    open (my $fh, $gtf_file) or die "Error, cannot open file $gtf_file";
    while (<$fh>) {
        chomp;
        my @x = split(/\t/);
        my $scaffold_id = $x[0];
        my $type = $x[2];
        
        unless ($type eq 'exon') { next; }
        

        my $interval_tree = $scaffold_itree_to_exon_structs_href->{$scaffold_id};
        unless (ref $interval_tree) {
            $interval_tree = $scaffold_itree_to_exon_structs_href->{$scaffold_id} = Set::IntervalTree->new;
        }
        
        my $info = $x[8];
        $info =~ /FI_gene_label \"([^\"]+)\"/ or die "Error, cannot parse FI_gene_label from $info";
        my $gene_id = $1;

        my ($lend, $rend) = ($x[3], $x[4]);
        
        
        # get original coordinate mapping info
        $info =~ /orig_coord_info \"([^,]+),(\d+),(\d+),([+-])\"/ or die "Error, cannot parse original coordinate info from $info";
        my $orig_chr = $1;
        my $orig_lend = $2;
        my $orig_rend = $3;
        my $orig_orient = $4;
        
        my ($orig_end5, $orig_end3) = ($orig_orient eq '+') ? ($orig_lend, $orig_rend) : ($orig_rend, $orig_lend);
        
        my $exon_struct = { scaffold_id => $scaffold_id,

                            gene_id => $gene_id,
                            
                            lend => $lend,
                            rend => $rend,
                
                            orig_chr => $orig_chr,
                            orig_end5 => $orig_end5,
                            orig_end3 => $orig_end3,
                            orig_orient => $orig_orient,

        };

        if ($lend != $rend) {
            $interval_tree->insert($exon_struct, $lend, $rend);
        }
        
        push (@{$scaff_to_gene_to_coords{$scaffold_id}->{$gene_id}}, $lend, $rend);
        
    }
    close $fh;

    
    my %scaffold_to_gene_structs;

    foreach my $scaffold (keys %scaff_to_gene_to_coords) {
        my @genes = keys %{$scaff_to_gene_to_coords{$scaffold}};
    
        my @gene_structs;
    
        foreach my $gene (@genes) {
            my @coords = sort {$a<=>$b} @{$scaff_to_gene_to_coords{$scaffold}->{$gene}};
            my $lend = shift @coords;
            my $rend = pop @coords;
            push (@{$scaffold_to_gene_structs{$scaffold}}, { gene_id => $gene,
                                                             lend => $lend,
                                                             rend => $rend,
                  });
        }
        
    }
    
    return(%scaffold_to_gene_structs);
}

####
sub divide_junction_read_at_breakpoint {
    my ($genome_coords_aref, $read_coords_aref, $brkpt_A) = @_;

    my @left_read_segments;
    my @right_read_segments;

    while (@$genome_coords_aref) {
        my $g_coords_aref = shift @$genome_coords_aref;
        my $r_coords_aref = shift @$read_coords_aref;
        
        my ($genome_lend, $genome_rend) = @$g_coords_aref;
        my ($read_lend, $read_rend) = @$r_coords_aref;

        if ($genome_rend <= $brkpt_A) {
            push (@left_read_segments, $r_coords_aref);
        }
        else {
            push (@right_read_segments, $r_coords_aref);
        }

    }

    return(\@left_read_segments, \@right_read_segments);
    
}


####
sub alignment_has_excessive_soft_clipping {
    my ($cigar, $max_end_clip) = @_;
    
    ## check left soft clip
    if ($cigar =~ /^(\d+)[SH]/) {
        my $clip_len = $1;
        if ($clip_len > $MAX_END_CLIP) {
            return(1);
        }
    }
    
    ## check right soft clip
    if ($cigar =~ /(\d+)[SH]$/) {
        my $clip_len = $1;
        if ($clip_len > $MAX_END_CLIP) {
            return(1);
        }
    }


    return(0); #ok
}

####
sub has_min_anchor_seq_entropy {
    my ($read_seq_sref, $segment_aref, $min_entropy) = @_;

    my ($read_lend, $read_rend) = sort {$a<=>$b} @$segment_aref; # sorting is important as (-) strand alignment reads are ordered oppositely
    
    
    my $read_subseq = substr($$read_seq_sref, $read_lend - 1, $read_rend - $read_lend + 1);

    my $entropy = &SeqUtil::compute_entropy($read_subseq);
    
    
    if ($entropy >= $min_entropy) {
        #print "$read_subseq : $entropy : OK\n";
        return(1);
    }
    else {
        #print "$read_subseq : $entropy : FAIL\n";
        return(0);
    }
}


####
sub partition_gene_structs {
    my (@gene_structs) = @_;

    # should already be sorted by lend.

    my $left_gene_struct = shift @gene_structs;
    my @left_gene_structs = ($left_gene_struct);

    my $left_gene_id = $left_gene_struct->{gene_id};
    my $left_gene_symbol = $left_gene_id;
    $left_gene_symbol =~ s/\^.*$//;

    my @right_gene_structs;

    foreach my $gene_struct (@gene_structs) {
        my $gene_id = $gene_struct->{gene_id};
        my ($gene_sym, $rest) = split(/\^/, $gene_id);

        if ($gene_sym eq $left_gene_symbol) {
            push (@left_gene_structs, $gene_struct);
        }
        else {
            push (@right_gene_structs, $gene_struct);
        }
    }

    unless (@left_gene_structs && @right_gene_structs) {
        die "Error, couldn't partition gene structs across contig breakpoint: " . Dumper(\@gene_structs);
    }
    
    @left_gene_structs = sort {$a->{rend}<=>$b->{rend}} @left_gene_structs;
    @right_gene_structs = sort {$a->{lend}<=>$b->{lend}} @right_gene_structs;

    # want those that are most adjacent to the central fusion contig join
    my @gene_struct_tuple = ($left_gene_structs[$#left_gene_structs], $right_gene_structs[0]);

    return(@gene_struct_tuple);
}

                                
    

    
