#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use SAM_reader;
use SAM_entry;
use Carp;
use Data::Dumper;
use Set::IntervalTree;
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);


my $MIN_SMALL_ANCHOR = 12;
my $MAX_MISMATCHES = 2;
my $MIN_LARGE_ANCHOR = 25;
my $MAX_END_CLIP = 3;

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
#  --MAX_MISMATCHES <int>     default: $MAX_MISMATCHES
#  --MIN_SMALL_ANCHOR <int>   default: $MIN_SMALL_ANCHOR
#  --MIN_LARGE_ANCHOR <int>   default: $MIN_LARGE_ANCHOR
#  --MAX_END_CLIP <int>       default: $MAX_END_CLIP
#
##############################################################


__EOUSAGE__

    ;


my $gtf_file;
my $bam_file;
my $help_flag;

&GetOptions('help|h' => \$help_flag,
            'gtf_file=s' => \$gtf_file,
            'bam_file=s' => \$bam_file,

            'MAX_MISMATCHES=i' => \$MAX_MISMATCHES,
            'MAX_END_CLIP=i' => \$MAX_END_CLIP,
            'MIN_SMALL_ANCHOR=i' => \$MIN_SMALL_ANCHOR,
            'MIN_LARGE_ANCHOR=i' => \$MIN_LARGE_ANCHOR,
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
    

    my %fusion_junctions;
    my %fusion_large_anchors;
    
    my $counter = 0;
    ## find the reads that matter:
    my $sam_reader = new SAM_reader($bam_file);
    while (my $sam_entry = $sam_reader->get_next()) {
            
        $counter++;
        if ($counter % 10000 == 0) { 
            print STDERR "\r[$counter]   ";
        }

        ## examine number of mismatches in read alignment
        my $line = $sam_entry->get_original_line();
        if ($line =~ /NM:i:(\d+)/) {
            my $mismatch_count = $1;
            if ($mismatch_count > $MAX_MISMATCHES) {
                next;
            }
        }
        
        ## check end clipping of alignment
        my $cigar = $sam_entry->get_cigar_alignment();
        if ($cigar =~ /^(\d+)[SH]/) {
            my $clip_len = $1;
            if ($clip_len > $MAX_END_CLIP) {
                next;
            }
        }
        
        my $read_name = $sam_entry->reconstruct_full_read_name();

        my $scaffold = $sam_entry->get_scaffold_name();
        if ($scaffold ne $prev_scaff) {
            
            if (! exists $scaffold_to_gene_structs{$scaffold}) { next; } # STAR aligns to the whole genome + fusion contigs.

            my @gene_structs = @{$scaffold_to_gene_structs{$scaffold}};
            if (scalar @gene_structs != 2) {
                die "Error, didn't find only 2 genes in the gtf file: " . Dumper(\@gene_structs);
            }
            
            @gene_structs = sort {$a->{lend} <=> $b->{lend}} @gene_structs;
            
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
                $genes_matched{ $exon_hit->{exon_struct}->{gene_id} }++;
            }
            my $num_genes_matched = scalar(keys %genes_matched);
            print "Genes matched: $num_genes_matched\n" if $DEBUG;
            unless ($num_genes_matched > 1) { next; }
            
                        
            my $junction_coord_token = &get_junction_coord_token(@exon_hits);
            if ($junction_coord_token) {
                
                my ($gene_A, $brkpt_A, $gene_B, $brkpt_B, $splice_token) = split(/\t/, $junction_coord_token);
                my ($left_brkpt_length, $right_brkpt_length) = &divide_junction_read_at_breakpoint($genome_coords_aref, $read_coords_aref, $brkpt_A);
                
                print STDERR "LEFT_read_seg: $left_brkpt_length, RIGHT_read_seg: $right_brkpt_length\n" if $DEBUG;

                if ($left_brkpt_length >= $MIN_SMALL_ANCHOR
                    &&
                    $right_brkpt_length >= $MIN_SMALL_ANCHOR) {
                    # calling it a fusion read.
                    
                    push (@{$fusion_junctions{$junction_coord_token}}, $read_name);
                    print $sam_entry->get_original_line() . "\n";
                
                    if ($left_brkpt_length >= $MIN_LARGE_ANCHOR) {
                        $fusion_large_anchors{$junction_coord_token}->{LEFT_LARGE_ANCHOR}++;
                    }
                    if ($right_brkpt_length >= $MIN_LARGE_ANCHOR) {
                        $fusion_large_anchors{$junction_coord_token}->{RIGHT_LARGE_ANCHOR}++;
                    }
                    
                }
                
            }
        }
        
    }
    
    {
        
        print STDERR "-writing fusion junction support info.\n";
        
        ## output the junction info summary:
        my $junction_coords_info_file = "$bam_file.fusion_junction_info";
        open (my $ofh, ">$junction_coords_info_file") or die "Error, cannot write to file: $junction_coords_info_file";
        
        my @fusion_j = reverse sort { $#{$fusion_junctions{$a}} <=> $#{$fusion_junctions{$b}} } keys %fusion_junctions;
        foreach my $fusion (@fusion_j) {
            my @reads = @{$fusion_junctions{$fusion}};
            my $num_fusion_reads = scalar @reads;

            my $LARGE_ANCHOR_BOTH_ENDS = ($fusion_large_anchors{$fusion}->{LEFT_LARGE_ANCHOR}
                                          &&
                                          $fusion_large_anchors{$fusion}->{RIGHT_LARGE_ANCHOR})
                ? "YES"
                : "NO";
            
            
            print $ofh "$fusion\t$num_fusion_reads\t$LARGE_ANCHOR_BOTH_ENDS\t" . join(",", @reads) . "\n";
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
        
        for (my $j = $i+1; $j <= $#exon_hits; $j++) {
            
            my $exon_hit_j = $exon_hits[$j];
            my $exon_struct_j = $exon_hit_j->{exon_struct};
            

            unless ($exon_struct_i->{rend} < $exon_struct_j->{lend}) { next; }
            
            ## see if different genes and properly split alignment
            if ($exon_struct_i->{gene_id} ne $exon_struct_j->{gene_id}
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
        $info =~ /gene_name \"([^\"]+)\"/ or die "Error, cannot parse gene_name from $info";
        my $gene_id = $1;

        my ($lend, $rend) = ($x[3], $x[4]);
        
        
        # get original coordinate mapping info
        $info =~ /orig_coord_info \"(\w+),(\d+),(\d+),([+-])\"/ or die "Error, cannot parse original coordinate info from $info";
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

    my $left_read_length = 0;
    my $right_read_length = 0;

    while (@$genome_coords_aref) {
        my $g_coords_aref = shift @$genome_coords_aref;
        my $r_coords_aref = shift @$read_coords_aref;
        
        my ($genome_lend, $genome_rend) = @$g_coords_aref;
        my ($read_lend, $read_rend) = @$r_coords_aref;

        my $read_segment_length = abs($read_rend - $read_lend) + 1;
        
        if ($genome_rend <= $brkpt_A) {
            $left_read_length += $read_segment_length;
        }
        else {
            $right_read_length += $read_segment_length;
        }

    }

    return($left_read_length, $right_read_length);
    
}
