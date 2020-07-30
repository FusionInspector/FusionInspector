#!/usr/bin/env perl

use strict;
use warnings;
use Carp;

use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use SAM_reader;
use SAM_entry;
use DelimParser;
use SeqUtil;
use TiedHash;
use JSON::XS;
use Overlap_piler;
use Data::Dumper;
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);


my $gtf_file;
my $bam_file;
my $junction_info_file;
my $genome_lib_dir;


my $MAX_END_CLIP = 10;

my $MIN_ALIGN_PER_ID = 96;
my $MIN_SEQ_ENTROPY = 1.2;

my $FUZZ = 5; # small fuzzy alignment bounds for spanning frags around breakpoint

my $usage = <<__EOUSAGE__;

###############################################################
#
# Required:
#
#  --gtf_file <string>         genePairContig.gtf
#  --bam <string>              read_alignments.bam
#  --junction_info <string>    bam.fusion_junction_info
#  --genome_lib_dir <string>   genome_lib_dir
#
# Optional:
#
#  --MIN_ALIGN_PER_ID <int>     default: $MIN_ALIGN_PER_ID
#  --MAX_END_CLIP <int>         default: $MAX_END_CLIP
#  --MIN_SEQ_ENTROPY <float>    default: $MIN_SEQ_ENTROPY
#
#  --debug|d
#
##############################################################



__EOUSAGE__

    ;

my $help_flag;
my $DEBUG = 0;


&GetOptions('help|h' => \$help_flag,
            
            'gtf_file=s' => \$gtf_file,
            'bam=s' => \$bam_file,
            'junction_info=s' => \$junction_info_file,
            'genome_lib_dir=s' => \$genome_lib_dir,

            'MIN_ALIGN_PER_ID=i' => \$MIN_ALIGN_PER_ID,
            'MAX_END_CLIP=i' => \$MAX_END_CLIP,
            'MIN_SEQ_ENTROPY=f' => \$MIN_SEQ_ENTROPY,

            'debug|d' => \$DEBUG,
    );

if ($help_flag) {
    die $usage;
}

unless ($gtf_file && $bam_file && $junction_info_file && $genome_lib_dir) {
    die $usage;
}


my %exon_bounds;
my %orig_coord_info;
my %scaffold_to_gene_structs;
my %scaffold_to_gene_breaks;
my %spanning_only_info;

my %junction_reads_ignore;
my %fusion_junctions;
my %fusion_breakpoint_info;

my %fusion_to_spanning_reads;
my %fusion_to_contrary_support;
my %spanning_read_want;



my $BLAST_ALIGNS_IDX;
my $blast_aligns_idx_file = "$genome_lib_dir/trans.blast.align_coords.align_coords.dbm";
if (-s $blast_aligns_idx_file) {
    $BLAST_ALIGNS_IDX = new TiedHash( { use => $blast_aligns_idx_file } );
}
else {
    die "Error, cannot locate blast idx file: $blast_aligns_idx_file";
}

my $JSON_DECODER = JSON::XS->new();



main: {


    ##############################################################
    ## Get the reference gene coordinates on each fusion-scaffold

    %scaffold_to_gene_structs = &parse_gtf_file($gtf_file, \%exon_bounds, \%orig_coord_info);

    {
        foreach my $scaffold (keys %scaffold_to_gene_structs) {

            my @gene_structs = @{$scaffold_to_gene_structs{$scaffold}};
            
            @gene_structs = sort {$a->{lend} <=> $b->{lend}} @gene_structs;

            if (scalar @gene_structs > 2) {
                @gene_structs = &partition_gene_structs(@gene_structs);
            }
                        
            
            my $left_gene = $gene_structs[0];
            my $right_gene = $gene_structs[1];
            
            my ($geneA, $geneB) = ($left_gene->{gene_id}, $right_gene->{gene_id});
            
            my $gene_bound_left = $left_gene->{rend};
            my $gene_bound_right = $right_gene->{lend};
            
            
            if ($gene_bound_left > $gene_bound_right) {
                die "Error, bounds out of order: $gene_bound_left ! <  $gene_bound_right";
            }
            $scaffold_to_gene_breaks{$scaffold} = [$gene_bound_left, $gene_bound_right];

            # for later, just in case we need to report info for those fusions lacking defined breakpoints
            $spanning_only_info{$scaffold} = join("\t", 
                                                  $geneA, $gene_bound_left, $orig_coord_info{$scaffold}->{$gene_bound_left}, 
                                                  $geneB, $gene_bound_right, $orig_coord_info{$scaffold}->{$gene_bound_right}, 
                                                  "NO_JUNCTION_READS_IDENTIFIED");
            
        }
    }
    

    
    #####################################
    ## Capture the fusion breakpoint info
    
    {
        open (my $fh, "$junction_info_file") or die "error, cannot open file $junction_info_file";
        my $tabreader = new DelimParser::Reader($fh, "\t");
        
        
        while (my $row = $tabreader->get_row()) {
            
            my ($geneA, $coordA, $orig_coordA, 
                $geneB, $coordB, $orig_coordB, 
                $splice_info, $fusion_read_count, 
                $large_anchor_support,
                $fusion_read_list) = ($row->{LeftGene},
                                      $row->{LeftLocalBreakpoint},
                                      $row->{LeftBreakpoint},
                                      $row->{RightGene},
                                      $row->{RightLocalBreakpoint},
                                      $row->{RightBreakpoint},
                                      $row->{SpliceType},
                                      $row->{JunctionReadCount},
                                      $row->{LargeAnchorSupport},
                                      $row->{JunctionReads});
            
            foreach my $fusion_read (split(/,/, $fusion_read_list)) {
                $fusion_read =~ s/\/[12]$//; #want core read name. 
                $junction_reads_ignore{$fusion_read}++;
            }

            
            my $gene_symA = $geneA;
            $gene_symA =~ s/\^.*$//;

            my $gene_symB = $geneB;
            $gene_symB =~ s/\^.*$//;
            
            my $fusion_name = join("--", $gene_symA, $gene_symB); # must match fusion scaffold name
            push (@{$fusion_junctions{$fusion_name}}, "$coordA-$coordB");
            
            my $breakpoint = "$coordA-$coordB";
            
            $fusion_breakpoint_info{"$fusion_name|$breakpoint"} = join("\t", $geneA, $coordA, $orig_coordA,
                                                                       $geneB, $coordB, $orig_coordB, $splice_info);
            
        }
        close $fh;
    

        if ($DEBUG) {
            print STDERR "Fusion breakpoint info: " . Dumper(\%fusion_breakpoint_info);
        }
    }


    ####################################################
    ## for each paired read, get the bounds of that read
    
       
    {

        #########################################################################
        ##  Generate the spanning fragment and contrary read support info summary
        

        ## output the spanning read info
        my $spanning_read_info_file = "$bam_file.fusion_spanning_info";
        open (my $ofh, ">$spanning_read_info_file") or die "Error, cannot write to $spanning_read_info_file";

        print STDERR "-outputting the spanning read info: $spanning_read_info_file.\n";
                
        my @fields = qw(LeftGene LeftLocalBreakpoint LeftBreakpoint
                        RightGene RightLocalBreakpoint RightBreakpoint
                        SpliceType 
                        SpanningFragCount SpanningFrags
                        NumCounterFusionLeft CounterFusionLeftReads
                        NumCounterFusionRight CounterFusionRightReads);

        my $tab_writer = new DelimParser::Writer($ofh, "\t", \@fields);



        my %scaffold_read_pair_to_read_bounds;
        my %core_counter;
        ## find the reads that matter:
        my $scaffold;
        my $prev_scaffold = "";
 
        my %read_alignment_counter = &count_read_alignments_among_fusion_contigs($bam_file);
        

        my %filtered_read_reason_counter;
        
        my $counter = 0;
        my $sam_reader = new SAM_reader($bam_file);
        while (my $sam_entry = $sam_reader->get_next()) {
            $counter++;
            print STDERR "\r[$counter]   " if $counter % 1000 == 0;
            
            if ($sam_entry->is_duplicate()) {
                next;
            }
            
            $scaffold = $sam_entry->get_scaffold_name();
            
            unless (exists $scaffold_to_gene_breaks{$scaffold}) { next; } # StarFI includes the whole genome, not just the fusion scaffs
            
            if ($scaffold ne $prev_scaffold) {
                if ($DEBUG) {
                    print STDERR "scaffold read pair to read bounds: " . Dumper(\%scaffold_read_pair_to_read_bounds);
                }
                
                if (%scaffold_read_pair_to_read_bounds) {
                    
                    &capture_spanning_frags($prev_scaffold, 
                                            \%scaffold_read_pair_to_read_bounds,
                                            \%core_counter,
                                            $tab_writer);
                    
                }
                %scaffold_read_pair_to_read_bounds = (); # reinit
                %core_counter = (); # reinit
                $prev_scaffold = $scaffold;
            }
            
            my $qual_val = $sam_entry->get_mapping_quality();
            
            my $read_name = $sam_entry->get_read_name();

            ## examine number of mismatches in read alignment
            my $line = $sam_entry->get_original_line();
            my $mismatch_count = 0;
            if ($line =~ /NM:i:(\d+)/i) {
                $mismatch_count = $1;
            }

            my $read_group;
            if ($line =~ /RG:Z:(\S+)/) {
                $read_group = $1;
            }
            
            my $alignment_length = $sam_entry->get_alignment_length();
            unless ($alignment_length) {
                if ($DEBUG) {
                    print STDERR "-skipping $read_name, no alignment length.\n";
                }
                next;
            }
            my $per_id = ($alignment_length - $mismatch_count) / $alignment_length * 100;
            if ($per_id < $MIN_ALIGN_PER_ID) {
                if ($DEBUG) {
                    print STDERR "-skipping $read_name, per_id $per_id < $MIN_ALIGN_PER_ID required.\n";
                }
                next;
            }

            $line =~ /NH:i:(\d+)/ or die "Error, cannot extract hit count (NH:i:) from $line";
            my $hit_count = $1;
                        
            ## check end clipping of alignment
            my $cigar = $sam_entry->get_cigar_alignment();
            if ($scaffold !~ /IGH/ && # //FIXME:  IGH here is a hack... should have more prinicpled ways of dealing with except cases.
                &alignment_has_excessive_soft_clipping($cigar)) {
                if ($DEBUG) {
                    print STDERR "-skipping $read_name, excessive softclipping: $cigar.\n";
                }
                next;
            }
                        
            my $scaffold_pos = $sam_entry->get_scaffold_position();

            my $mate_scaffold_name = $sam_entry->get_mate_scaffold_name();            
            my $mate_scaffold_pos = $sam_entry->get_mate_scaffold_position();
            
            unless ($mate_scaffold_name eq $scaffold || $mate_scaffold_name eq "=") { next; }


            my $core_read_name = $sam_entry->get_core_read_name();            

            if ($junction_reads_ignore{$core_read_name}) { 
                # junction reads cannot be used as spanning frags.
                if ($DEBUG) {
                    print STDERR "-skipping $read_name, prev identified as a breakpoint (junction) read.\n";
                }
                next; 
            } 


            # check if alignments begin in their respective fusion gene areas:
            my ($scaff_gene_left_rend, $scaff_gene_right_lend) = @{$scaffold_to_gene_breaks{$scaffold}};
            
            my ($span_lend, $span_rend) = sort {$a<=>$b} $sam_entry->get_genome_span();
            
            # be sure that each read on its own is entirely encapsulated within a single gene region (not crossing the bounds)
            unless ($span_rend <= $scaff_gene_left_rend || $span_lend >= $scaff_gene_right_lend) { 
                
                if ($DEBUG) {
                    print STDERR "-skipping $read_name with span [$span_lend-$span_rend], not restricted to one side of fusion inter-region [$scaff_gene_left_rend-$scaff_gene_right_lend]\n";
                }
                
                next; 
            }
            
            my $alignment_side = ($span_rend <= $scaff_gene_left_rend) ? "LEFT" : "RIGHT";
            
            my $strand = $sam_entry->get_query_strand();
            
            my $token = join("$;", $read_name, $scaffold);

            
            my $read_seq = $sam_entry->get_sequence();
            my $entropy = &SeqUtil::compute_entropy($read_seq);
            unless ($entropy >= $MIN_SEQ_ENTROPY) { 
                if ($DEBUG) {
                    print STDERR "-skipping $read_name, entropy $entropy < $MIN_SEQ_ENTROPY required.\n" if $DEBUG;
                    $filtered_read_reason_counter{"low entroy"} += 1;
                }
                next; 
            }
            
            
            my ($genome_coords_aref, $read_coords_aref) = $sam_entry->get_alignment_coords();
            my @align_segment_overlap_pairs = &get_overlapping_exons($genome_coords_aref, $exon_bounds{$scaffold});
            if (@align_segment_overlap_pairs) {
                
                ## check if in seq-similar regions between gene pairs.
                ##
                
                if (&exceedingly_overlaps_homologous_segment($scaffold, $alignment_side, $genome_coords_aref, \@align_segment_overlap_pairs, $orig_coord_info{$scaffold})) {
                    print STDERR "-skipping $read_name, aligns to seq-similar contig region between gene pairs\n" if $DEBUG;
                    $filtered_read_reason_counter{"seq similar region alignment"} += 1;
                    next;
                }
                
            }
            else {
                
                #print STDERR "No exon overlap: " . Dumper($genome_coords_aref) . Dumper($exon_bounds{$scaffold});
                $filtered_read_reason_counter{"lacks exon overlap"} += 1;
                next; 
            } # only examine exon-overlapping entries
            
            

            my $full_read_name = $sam_entry->reconstruct_full_read_name();
            if ($full_read_name =~ /^(\S+)\/([12])$/) {
                my ($core, $pair_end) = ($1, $2);
                $core_counter{"$scaffold|$core"}++;  # track how many alignments we have for this rnaseq fragment
                
                $scaffold_read_pair_to_read_bounds{$scaffold}->{$core}->[$pair_end-1] = { span_lend => $span_lend, 
                                                                                          span_rend => $span_rend, 
                                                                                          strand => $strand,
                                                                                          qual => $qual_val,
                                                                                          full_read_name => $full_read_name,
                                                                                          NH => $hit_count,
                                                                                          fusion_scaff_hit_count => $read_alignment_counter{$full_read_name},
                                                                                          read_group => $read_group,
                };
                
                
            }
        }
        
        if (%scaffold_read_pair_to_read_bounds) {
            &capture_spanning_frags($scaffold, \%scaffold_read_pair_to_read_bounds, \%core_counter, $tab_writer);
        }
        
        close $ofh; # donw writing fusion report.

        print STDERR "-filtered reads reasons: " . Dumper(\%filtered_read_reason_counter);
        
    }
    

    #################################################
    # output the spanning reads we want in SAM format
    
    if (%spanning_read_want) {

        my $num_spanning_reads_want = scalar(keys %spanning_read_want);
        print STDERR "-retrieving read alignments for $num_spanning_reads_want spanning frags.\n";
        
        my %missing = %spanning_read_want;
        
        my $sam_reader = new SAM_reader($bam_file);
        while (my $sam_entry = $sam_reader->get_next()) {
                        
            my $scaffold = $sam_entry->get_scaffold_name();
            my $core_read_name = $sam_entry->get_core_read_name();
            if ($spanning_read_want{"$scaffold|$core_read_name"}) {
                print $sam_entry->get_original_line() . "\n";

                if (exists $missing{"$scaffold|$core_read_name"}) {
                    delete $missing{"$scaffold|$core_read_name"};
                }
            }
        }

        if (%missing) {
            confess "Error, didn't extract the following spanning frags from the bam file ($bam_file): " . Dumper(\%missing);
        }
    }
    

    exit(0);
}


####
sub get_overlapping_exons {
    my ($genome_coords_aref, $exon_bounds_href) = @_;

    my @overlapping_segment_pairs;

    foreach my $coordset (@$genome_coords_aref) {
        my ($lend, $rend) = @$coordset;
        foreach my $exon_coordset (keys %$exon_bounds_href) {
            my ($e_lend, $e_rend) = split(/-/, $exon_coordset);
            if ($e_lend < $rend && $e_rend > $lend) {
                my $exon_coordset_aref = [$e_lend, $e_rend];
                push (@overlapping_segment_pairs, [$coordset, $exon_coordset_aref]);
            }
        }
    }
    
    return(@overlapping_segment_pairs);
}



####
sub parse_gtf_file {
    my ($gtf_file, $exon_bounds_href, $orig_coord_info_href) = @_;

    my %scaff_to_gene_to_coords;

    open (my $fh, $gtf_file) or die "Error, cannot open file $gtf_file";
    while (<$fh>) {
        chomp;
        my @x = split(/\t/);
        my $scaffold_id = $x[0];
        my $type = $x[2];
        
        unless ($type eq 'exon') { next; }
        
        my $info = $x[8];
        $info =~ /FI_gene_label \"([^\"]+)\"/ or die "Error, cannot parse FI_gene_label from $info";
        my $gene_id = $1;

        my ($lend, $rend) = ($x[3], $x[4]);
        push (@{$scaff_to_gene_to_coords{$scaffold_id}->{$gene_id}}, $lend, $rend);
        $exon_bounds_href->{$scaffold_id}->{"$lend-$rend"} = 1;
        
        #  orig_coord_info "chr7,34697897,34698171,+";
        $info =~ /orig_coord_info \"([^\"]+)\"/ or die "Error, cannot parse orig_coord_info from $info";
        my $orig_coord_info = $1;
        my ($orig_chr, $orig_lend, $orig_rend, $orig_orient) = split(/,/, $orig_coord_info);
        $orig_coord_info_href->{$scaffold_id}->{$lend} = join(":", $orig_chr, $orig_lend, $orig_orient);
        $orig_coord_info_href->{$scaffold_id}->{$rend} = join(":", $orig_chr, $orig_rend, $orig_orient);
        
        

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


####
sub capture_spanning_frags {
    my ($scaffold, $scaffold_read_pair_to_read_bounds_href, $core_counter_href, $tab_writer) = @_;
    
    print STDERR "-fusion SPANNING read extraction for scaff: $scaffold\n";
    
    my %scaffold_read_pair_to_read_bounds = %$scaffold_read_pair_to_read_bounds_href;
    my %core_counter = %$core_counter_href;
    
    ##########################################
    # determine which reads are spanning reads
    
    my ($gene_bound_left, $gene_bound_right) = @{$scaffold_to_gene_breaks{$scaffold}};
                        
    if ($gene_bound_left > $gene_bound_right) { 
        die "Error, gene bounds out of range for $scaffold: $gene_bound_left - $gene_bound_right "; 
    }
    
    
    foreach my $fragment (keys %{$scaffold_read_pair_to_read_bounds{$scaffold}}) {
        
        if ($core_counter{"$scaffold|$fragment"} != 2) { next; } # ignore those fragments that have multiply-mapping reads to this contig.
        
        my @pair_coords = grep { defined $_ } @{$scaffold_read_pair_to_read_bounds{$scaffold}->{$fragment}};
        
        ## data structure format: $scaffold_read_pair_to_read_bounds{$scaffold}->{$core}->[$pair_end-1] = [$span_lend, $span_rend, $strand];
        
        if (scalar @pair_coords > 1) {
            # need both paired ends
            
            @pair_coords = sort {$a->{span_lend} <=> $b->{span_lend}} @pair_coords;
            
            my $left_read_lend = $pair_coords[0]->{span_lend};
            my $left_read_rend = $pair_coords[0]->{span_rend};
            
            my $right_read_lend = $pair_coords[1]->{span_lend};
            my $right_read_rend = $pair_coords[1]->{span_rend};
            
            my $left_read_orient = $pair_coords[0]->{strand};
            my $right_read_orient = $pair_coords[1]->{strand};
            
            my $read_group = $pair_coords[0]->{read_group};
            
            unless ($left_read_orient eq '+' && $right_read_orient eq '-') { 
                # not proper pairs after all!!
                print STDERR "-skipping pair $scaffold|$fragment, discordantly aligned\n" if $DEBUG;
                next; 
            }  
            
            #####################################################################
            ## assign fragment as fusion support based on breakpoint coordinates.
            
            
            my $is_fusion_spanning_fragment_flag = 0;
            if ($left_read_rend <= $gene_bound_left && $right_read_lend >= $gene_bound_right
                
                # ensure ok quality
                # && $pair_coords[0]->{qual} > 0 && $pair_coords[1]->{qual} > 0
                
                # ensure single hit
                && $pair_coords[0]->{NH} == $pair_coords[0]->{fusion_scaff_hit_count} 
                && $pair_coords[1]->{NH} == $pair_coords[1]->{fusion_scaff_hit_count} # yes, must be uniquely supporting the fusion here!
                ) 
            {
                $is_fusion_spanning_fragment_flag = 1;
                $spanning_read_want{"$scaffold|$fragment"}++; # capture for SAM-retreival next.
                #print STDERR "-want $scaffold|$fragment in sam\n";
            }
            else {
                if ($DEBUG) {
                    print STDERR "-fragment: $scaffold|$fragment not flagged as fusion spanning. "
                        . " frag coords: $left_read_lend-$left_read_rend:$left_read_orient -- $right_read_lend-$right_read_rend:$right_read_orient gene_bounds: $gene_bound_left,$gene_bound_right\n" . Dumper(\@pair_coords);
                }
            }
            
            
            ##########
            ## encode the read group into the fragment name:
            
            if ($read_group) {
                $fragment = "&" . $read_group . "@" . $fragment;
            }
            
            
            #################
            ## assign spanning frags to the specific breakpoints
            
            my $assigned_to_breakpoint_flag = 0;
            
            my $candidate_fusion_breakpoints_aref = $fusion_junctions{$scaffold};
            if (ref $candidate_fusion_breakpoints_aref) {
                
                print STDERR "Candidate fusion breakpoints for scaffold: $scaffold: " . Dumper($candidate_fusion_breakpoints_aref) if $DEBUG;
                
                foreach my $fusion_breakpoint (@{$candidate_fusion_breakpoints_aref}) {
                    my ($break_lend, $break_rend) = split(/-/, $fusion_breakpoint);
                    print STDERR "** Breakpoint: $break_lend, $break_rend\n" if $DEBUG;
                    
                    print STDERR "$fragment\tr1: $left_read_lend-$left_read_rend   r2: $right_read_lend-$right_read_rend  brk: $fusion_breakpoint\n" if $DEBUG;
                    
                    if ($is_fusion_spanning_fragment_flag) {
                        
                        if ($left_read_rend <= $break_lend + $FUZZ && $break_rend - $FUZZ < $right_read_lend) {
                            
                            # <=======>                                                    <=======>   # reads
                            #           |------------------------------------------------|   # breakpoints on scaffold
                            
                            
                            # junction-specific spanning fragment support assignment
                            
                            # must meet the more restrictive criteria wrt qual and NH
                            push (@{$fusion_to_spanning_reads{"$scaffold|$fusion_breakpoint"}}, $fragment);
                            $assigned_to_breakpoint_flag = 1;
                            
                            #print STDERR "\t-capturing spanning breakpoint frag: $fragment\n";
                        }
                        
                    }
                    else {
                        ## Not a fusion-spanning fragment.
                        ## See if it's fusion-countering evidence
                        
                        
                        if ($left_read_lend < $break_lend 
                            && $break_lend < $right_read_rend
                            && $right_read_rend < $gene_bound_right) {
                            
                            # <==---?------?----==>   # reads
                            #           |------------------------------------------------|   # breakpoints on scaffold
                            
                            
                            
                            ## contrary support at left junction
                            push (@{$fusion_to_contrary_support{"$scaffold|$fusion_breakpoint"}->{left}}, $fragment);
                        }
                        elsif ($left_read_lend < $break_rend 
                               && $break_rend < $right_read_rend
                               && $left_read_lend > $gene_bound_left) {
                            
                            
                            #                                                  <==-----?----?----===>   # reads
                            #           |------------------------------------------------|   # breakpoints on scaffold
                            
                            
                            
                            ## contrary support at right junction
                            push (@{$fusion_to_contrary_support{"$scaffold|$fusion_breakpoint"}->{right}}, $fragment);
                        }
                        
                    }
                } # end of foreach breakpoing candidate
            } # end of if have candidate breakpoints
            
            if ($is_fusion_spanning_fragment_flag && ! $assigned_to_breakpoint_flag) {
                # still capture it even though there's no junction read to assign it to.
                
                my $fuzzy_breakpoint = join("-", $gene_bound_left, $gene_bound_right);
                push (@{$fusion_to_spanning_reads{"$scaffold|$fuzzy_breakpoint"}}, $fragment);
                
            }
            
        } # endif have read pair coords
    } # end of foreach fragment
    
    
    
    ## output fusion records for that scaffold:
    foreach my $fusion_n_breakpoint (sort keys %fusion_to_spanning_reads) {
            my ($fusion_name, $breakpoint) = split(/\|/, $fusion_n_breakpoint);
            my ($geneA, $geneB) = split(/--/, $fusion_name);

            my $fusion_info = $fusion_breakpoint_info{$fusion_n_breakpoint};
            
            unless ($fusion_info) {
                # spanning frags, no junction breakpoint
                $fusion_info = $spanning_only_info{$fusion_name};
                unless ($fusion_info) {
                    confess "Error, no fusion info for [$fusion_name] ";
                }
            }
            
            my ($coordA, $coordB) = split(/-/, $breakpoint);
            
            my @spanning_reads = @{$fusion_to_spanning_reads{$fusion_n_breakpoint}};
            
            my $num_spanning = scalar(@spanning_reads);
            
            ## examine contrary support
            my @contrary_left_support;
            if (my $support_aref = $fusion_to_contrary_support{$fusion_n_breakpoint}->{left}) {
                @contrary_left_support = @$support_aref;
            }
            my @contrary_right_support;
            if (my $support_aref = $fusion_to_contrary_support{$fusion_n_breakpoint}->{right}) {
                @contrary_right_support = @$support_aref;
            }
            
            
            my $num_left_contrary_support = scalar(@contrary_left_support);
            my $num_right_contrary_support = scalar(@contrary_right_support);
            
            my $contrary_left_support_txt = join(",", @contrary_left_support) || ".";
            my $contrary_right_support_txt = join(",", @contrary_right_support) || ".";

            
            my ($LeftGene, $LeftLocalBreakpoint, $LeftBreakpoint,
                $RightGene, $RightLocalBreakpoint, $RightBreakpoint,
                $SpliceType) = split(/\t/, $fusion_info);

            
            $tab_writer->write_row( { LeftGene => $LeftGene,
                                      LeftLocalBreakpoint => $LeftLocalBreakpoint,
                                      LeftBreakpoint => $LeftBreakpoint,
                                      RightGene => $RightGene,
                                      RightLocalBreakpoint => $RightLocalBreakpoint,
                                      RightBreakpoint => $RightBreakpoint,
                                      SpliceType => $SpliceType,
                                      SpanningFragCount => $num_spanning,
                                      SpanningFrags => join(",", @spanning_reads),
                                      NumCounterFusionLeft => $num_left_contrary_support,
                                      CounterFusionLeftReads => $contrary_left_support_txt,
                                      NumCounterFusionRight => $num_right_contrary_support,
                                      CounterFusionRightReads => $contrary_right_support_txt,
                                    } );
            
        }



    return;
}
    

####
sub count_read_alignments_among_fusion_contigs {
    my ($bam_file) = @_;

    my %alignment_counter;

    my $sam_reader = new SAM_reader($bam_file);
    while (my $sam_entry = $sam_reader->get_next()) {
        
        # ensure on fusion contig.
        my $contig = $sam_entry->get_scaffold_name();
        unless ($contig =~ /\-\-/) {
            next;
        }
        
        my $full_read_name = $sam_entry->reconstruct_full_read_name();
        
        $alignment_counter{$full_read_name} += 1;
    }

    return(%alignment_counter);
}
    
        
####
sub exceedingly_overlaps_homologous_segment {
    my ($scaffold, $alignment_side, $read_genome_coords_aref, $align_segment_overlap_pairs_aref, $original_genome_coord_mapping_href) = @_;
    
    my $blast_pair_info = $BLAST_ALIGNS_IDX->get_value($scaffold);
    unless (defined $blast_pair_info) {
        return(0);
    }
    
    my $blast_align_info_struct = $JSON_DECODER->decode($blast_pair_info);
    
    my $seq_similar_genome_coords_aref = ($alignment_side eq "LEFT") ? $blast_align_info_struct->{'coords_A'} : $blast_align_info_struct->{'coords_B'};
    
    if ($DEBUG) {
        print STDERR "exceedingly_overlaps_homologous_segment: inputs:\n"
            . " scaffold: $scaffold\n"
            . " read_genome_coords_aref: " . Dumper($read_genome_coords_aref) . "\n"
            . " align_segment_overlap_pairs_aref: " . Dumper($align_segment_overlap_pairs_aref) . "\n"
            . " original_genome_coord_mapping_href: " . Dumper($original_genome_coord_mapping_href) . "\n";
    }
    
        
    my $read_segment_length = 0;
    foreach my $coordset (@$read_genome_coords_aref) {
        my ($lend, $rend) = sort {$a<=>$b} @$coordset;
        $read_segment_length += $rend - $lend + 1;
    }
    

    ## get read alignment coords in ref genome coord system
    
    my @read_genome_coords;
    
    foreach my $align_segment_overlap_pair (@$align_segment_overlap_pairs_aref) {
        
        my ($read_coordset, $exon_coordset) = @$align_segment_overlap_pair;
                
        ## remap the read coordinates to the genome reference:
        my ($new_read_genome_lend, $new_read_genome_rend) = sort {$a<=>$b} &convert_from_FI_contig_to_genome_coordinates($read_coordset,
                                                                                                                         $exon_coordset,
                                                                                                                         $original_genome_coord_mapping_href);
        
        push (@read_genome_coords, [$new_read_genome_lend, $new_read_genome_rend]);
        
    }
    
    ## check over overlapping regions w/ seq-similar segments

    my @overlapping_regions;
    foreach my $seq_similar_region (@$seq_similar_genome_coords_aref) {
        
        my ($seq_similar_genome_lend, $seq_similar_genome_rend) = sort {$a<=>$b} @$seq_similar_region;
        
        foreach my $read_genome_coordset (@read_genome_coords) {
            my ($read_lend, $read_rend) = @$read_genome_coordset;
            
            if ($seq_similar_genome_lend < $read_rend && $seq_similar_genome_rend > $read_lend) {
                my @coords = sort {$a<=>$b} ($read_lend, $read_rend, $seq_similar_genome_lend, $seq_similar_genome_rend);
                my ($overlap_lend, $overlap_rend) = ($coords[1], $coords[2]);
                push (@overlapping_regions, [$overlap_lend, $overlap_rend]);
            }
        }
        
    }
    
    print STDERR "-seqsimilar overlapping regions of read alignments: " . Dumper(\@overlapping_regions) if $DEBUG;

    if (@overlapping_regions) {
        
        my @collapsed_coords = &Overlap_piler::simple_coordsets_collapser(@overlapping_regions);
        
        print STDERR "-collapsed as: " . Dumper(\@collapsed_coords) if $DEBUG;
        
        my $seq_similar_region_len = 0;
        foreach my $collapsed_coordset (@collapsed_coords) {
            my ($lend, $rend) = @$collapsed_coordset;
            $seq_similar_region_len += $rend - $lend + 1;
        }
        
        if ($seq_similar_region_len / $read_segment_length > 0.5) {
            # majority of read segment is in seq similar region
            return(1);
        }
    }
    
    return(0);
        

}



####
sub convert_from_FI_contig_to_genome_coordinates {
    my ($read_coords_aref, $exon_coords_aref, $original_genome_coord_mapping_href) = @_;
            
    my ($read_FI_lend, $read_FI_rend) = sort {$a<=>$b} @$read_coords_aref;
    my ($exon_lend, $exon_rend) = @$exon_coords_aref;
    
    
    ##   genome_lend            genome_rend
    ##      |----------------------|  
    ##           |----------|
    ##       read_lend     read_rend
    
    ##     |-----|          |------|
    ##     genome_lend_delta   genome_rend_delta
    
    my $genome_lend_delta = $read_FI_lend - $exon_lend;
    my $genome_rend_delta = $exon_rend - $read_FI_rend;


    my $new_read_genome_lend;
    my $new_read_genome_rend;
    
    my $exon_end5 = (split(/:/, $original_genome_coord_mapping_href->{$exon_lend}))[1];
    my $exon_end3 = (split(/:/, $original_genome_coord_mapping_href->{$exon_rend}))[1];
    
    my $orig_orient = ($exon_end5 < $exon_end3) ? '+' : '-';
    
    ##   genome_end5   -->      genome_end3
    ##      |----------------------|  
    if ($orig_orient eq '+') {
            
        $new_read_genome_lend = $exon_end5 + $genome_lend_delta;
        $new_read_genome_rend = $exon_end3 - $genome_rend_delta;
    }
    
    
    ##   genome_end3  <--     genome_end5   ## revcomp orientation
    ##      |----------------------|  
    
    else {
        $new_read_genome_lend = $exon_end5 - $genome_lend_delta;
        $new_read_genome_rend = $exon_end3 + $genome_rend_delta;
    }
    
    return($new_read_genome_lend, $new_read_genome_rend);
    
}
