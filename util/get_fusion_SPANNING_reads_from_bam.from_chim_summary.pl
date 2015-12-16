#!/usr/bin/env perl

use strict;
use warnings;
use Carp;

use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use SAM_reader;
use SAM_entry;
use Data::Dumper;
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);


my $gtf_file;
my $bam_file;
my $junction_info_file;

my $MAX_END_CLIP = 10;
my $MAX_MISMATCHES = 2;

my $usage = <<__EOUSAGE__;

###############################################################
#
# Required:
#
#  --gtf_file <string>         genePairContig.gtf
#  --bam <string>              read_alignments.bam
#  --junction_info <string>    bam.fusion_junction_info
#
# Optional:
#
#  --MAX_MISMATCHES <int>     default: $MAX_MISMATCHES
#  --MAX_END_CLIP <int>       default: $MAX_END_CLIP
#
##############################################################

__EOUSAGE__

    ;

my $help_flag;

&GetOptions('help|h' => \$help_flag,
            
            'gtf_file=s' => \$gtf_file,
            'bam=s' => \$bam_file,
            'junction_info=s' => \$junction_info_file,
            
            'MAX_MISMATCHES=i' => \$MAX_MISMATCHES,
            'MAX_END_CLIP=i' => \$MAX_END_CLIP,
            
    );

if ($help_flag) {
    die $usage;
}

unless ($gtf_file && $bam_file) {
    die $usage;
}





main: {


    ##############################################################
    ## Get the reference gene coordinates on each fusion-scaffold
    

    my %exon_bounds;
    my %orig_coord_info;
    my %scaffold_to_gene_structs = &parse_gtf_file($gtf_file, \%exon_bounds, \%orig_coord_info);

    my %spanning_only_info;
    
    my %scaffold_to_gene_breaks;
    {
        foreach my $scaffold (keys %scaffold_to_gene_structs) {

            my ($geneA, $geneB) = split(/--/, $scaffold);

            my @gene_structs = @{$scaffold_to_gene_structs{$scaffold}};
            if (scalar @gene_structs != 2) {
                die "Error, didn't find only 2 genes in the gtf file: " . Dumper(\@gene_structs);
            }
            
            @gene_structs = sort {$a->{lend} <=> $b->{lend}} @gene_structs;
            
            my $left_gene = $gene_structs[0];
            my $right_gene = $gene_structs[1];
            
            my $gene_bound_left = $left_gene->{rend};
            my $gene_bound_right = $right_gene->{lend};
            
            
            if ($gene_bound_left > $gene_bound_right) {
                die "Error, bounds out of order: $gene_bound_left ! <  $gene_bound_right";
            }
            $scaffold_to_gene_breaks{$scaffold} = [$gene_bound_left, $gene_bound_right];

            $spanning_only_info{$scaffold} = join("\t", 
                                                  $geneA, $gene_bound_left, $orig_coord_info{$scaffold}->{$gene_bound_left}, 
                                                  $geneB, $gene_bound_right, $orig_coord_info{$scaffold}->{$gene_bound_right}, 
                                                  "NO_JUNCTION_READS_IDENTIFIED");
            
        }
    }
    
    #print STDERR "Scaffold to gene breaks: " . Dumper(\%scaffold_to_gene_breaks);
    

    
    my %junction_reads_ignore;
    my %fusion_junctions;
    my %fusion_breakpoint_info;

    
    {
        open (my $fh, "$junction_info_file") or die "error, cannot open file $junction_info_file";
        while (<$fh>) {
            chomp;
            my ($geneA, $coordA, $orig_coordA, $geneB, $coordB, $orig_coordB, $splice_info, $fusion_read_count, $fusion_read_list) = split(/\t/);

            foreach my $fusion_read (split(/,/, $fusion_read_list)) {
                $junction_reads_ignore{$fusion_read}++;
            }
            
            my $fusion_name = join("--", $geneA, $geneB);
            push (@{$fusion_junctions{$fusion_name}}, "$coordA-$coordB");
            
            my $breakpoint = "$coordA-$coordB";
            
            $fusion_breakpoint_info{"$fusion_name|$breakpoint"} = join("\t", $geneA, $coordA, $orig_coordA,
                                                                       $geneB, $coordB, $orig_coordB, $splice_info);
            
        }
        close $fh;
    }

    ## for each paired read, get the bounds of that read
    my %scaffold_read_pair_to_read_bounds;    
    my %core_counter;
    {
        ## find the reads that matter:
        my $counter = 0;
        my $sam_reader = new SAM_reader($bam_file);
        while (my $sam_entry = $sam_reader->get_next()) {
            $counter++;
            print STDERR "\r[$counter]   " if $counter % 1000 == 0;
            
            #print STDERR Dumper($sam_entry);

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
                        
            my $qual_val = $sam_entry->get_mapping_quality();
            #unless ($qual_val >= $MIN_QUALITY) { next; }
            

            my $scaffold = $sam_entry->get_scaffold_name();
            unless (exists $scaffold_to_gene_breaks{$scaffold}) { next; } # StarFI includes the whole genome, not just the fusion scaffs


            my $scaffold_pos = $sam_entry->get_scaffold_position();

            my $mate_scaffold_name = $sam_entry->get_mate_scaffold_name();
            my $mate_scaffold_pos = $sam_entry->get_mate_scaffold_position();
            
            unless ($mate_scaffold_name eq $scaffold || $mate_scaffold_name eq "=") { next; }

            # check if alignments begin in their respective fusion gene areas:
            my ($scaff_gene_left_rend, $scaff_gene_right_lend) = @{$scaffold_to_gene_breaks{$scaffold}};
            

            ## see that the fragment read pair together span both the genes
            my ($pos1, $pos2) = sort {$a<=>$b} ($scaffold_pos, $mate_scaffold_pos);
            
            unless ($pos1 < $scaff_gene_left_rend && $pos2 > $scaff_gene_right_lend) { next; }
            

            my ($span_lend, $span_rend) = sort {$a<=>$b} $sam_entry->get_genome_span();
            
            # be sure that each read on its own is entirely encapsulated within a single gene region (not crossing the bounds)
            unless ($span_rend <= $scaff_gene_left_rend || $span_lend >= $scaff_gene_right_lend) { next; }
            

            my $read_name = $sam_entry->get_read_name();
            
            my $strand = $sam_entry->get_query_strand();
            
            my $token = join("$;", $read_name, $scaffold);

            my $full_read_name = $sam_entry->reconstruct_full_read_name();            

            if ($junction_reads_ignore{$full_read_name}) { next; }
            
        
            my ($genome_coords_aref, $read_coords_aref) = $sam_entry->get_alignment_coords();
            unless (&overlaps_exon($genome_coords_aref, $exon_bounds{$scaffold}) ) { 
                #print STDERR "No exon overlap: " . Dumper($genome_coords_aref) . Dumper($exon_bounds{$scaffold});
                next; 
            } # only examine exon-overlapping entries
            
            
            if ($full_read_name =~ /^(\S+)\/([12])$/) {
                my ($core, $pair_end) = ($1, $2);
                $core_counter{"$scaffold|$core"}++;  # track how many alignments we have for this rnaseq fragment
                
                $scaffold_read_pair_to_read_bounds{$scaffold}->{$core}->[$pair_end-1] = [$span_lend, $span_rend, $strand];
                
            }
        }
    }
   
    my %fusion_to_spanning_reads;

    my %fusion_to_contrary_support;
    
    # determine which reads are spanning reads
    my %spanning_read_want;
    {
        foreach my $scaffold (keys %scaffold_to_gene_breaks) {
            my ($gene_bound_left, $gene_bound_right) = @{$scaffold_to_gene_breaks{$scaffold}};
            
            if ($gene_bound_left > $gene_bound_right) { 
                die "Error, gene bounds out of range for $scaffold: $gene_bound_left - $gene_bound_right "; 
            }
            
            
            foreach my $fragment (keys %{$scaffold_read_pair_to_read_bounds{$scaffold}}) {
                
                if ($core_counter{"$scaffold|$fragment"} != 2) { next; } # ignore those fragments that have multiply-mapping reads to this contig.
                my @pair_coords = grep { defined $_ } @{$scaffold_read_pair_to_read_bounds{$scaffold}->{$fragment}};
                if (scalar @pair_coords > 1) {
                    # need both paired ends
                    
                    @pair_coords = sort {$a->[0] <=> $b->[0]} @pair_coords;
                    
                    my $left_read_rend = $pair_coords[0]->[1];
                    my $right_read_lend = $pair_coords[1]->[0];
                    
                    my $left_read_orient = $pair_coords[0]->[2];
                    my $right_read_orient = $pair_coords[1]->[2];

                    unless ($left_read_orient eq '+' && $right_read_orient eq '-') { next; }  # not proper pairs after all!!

                    #####################################################################
                    ## assign fragment as fusion support based on breakpoint coordinates.


                    my $is_fusion_spanning_fragment_flag = 0;
                    if ($left_read_rend < $gene_bound_left && $right_read_lend > $gene_bound_right) {
                        $is_fusion_spanning_fragment_flag = 1;
                        $spanning_read_want{"$scaffold|$fragment"}++; # capture for SAM-retreival next.
                    }
    
                    my $candidate_fusion_breakpoints_aref = $fusion_junctions{$scaffold};
                    if (ref $candidate_fusion_breakpoints_aref) {
                                        
                        foreach my $fusion_breakpoint (@{$candidate_fusion_breakpoints_aref}) {
                            my ($break_lend, $break_rend) = split(/-/, $fusion_breakpoint);
                        
                            if ($is_fusion_spanning_fragment_flag) {
                                
                                ## assign all spanning frags to all junction reads.  Important for later filtering
                                push (@{$fusion_to_spanning_reads{"$scaffold|$fusion_breakpoint"}}, $fragment);
                                                              
                            }
                            else {
                                if ($left_read_rend < $break_lend && $break_lend < $right_read_lend
                                    && $right_read_lend < $gene_bound_right) {
                                    
                                    ## contrary support at left junction
                                    push (@{$fusion_to_contrary_support{"$scaffold|$fusion_breakpoint"}->{left}}, $fragment);
                                }
                                elsif ($left_read_rend < $break_rend && $break_rend < $right_read_lend
                                       && $left_read_rend > $gene_bound_left) {
                                    ## contrary support at right junction
                                    push (@{$fusion_to_contrary_support{"$scaffold|$fusion_breakpoint"}->{right}}, $fragment);
                                }
                        
                            }
                        }
                    }
                    elsif($is_fusion_spanning_fragment_flag) {
                        
                        # still capture it even though there's no junction read to assign it to.
                        
                        my $fuzzy_breakpoint = join("-", $gene_bound_left, $gene_bound_right);
                        push (@{$fusion_to_spanning_reads{"$scaffold|$fuzzy_breakpoint"}}, $fragment);
                        
                    }
                    
                    
                }
            }
        }
    }
    

    # output the spanning reads we want
    if (%spanning_read_want) {
        my $sam_reader = new SAM_reader($bam_file);
        while (my $sam_entry = $sam_reader->get_next()) {
            
            my $qual_val = $sam_entry->get_mapping_quality();
            #unless ($qual_val >= $MIN_QUALITY) { next; }
            
            
            my $scaffold = $sam_entry->get_scaffold_name();
            my $core_read_name = $sam_entry->get_core_read_name();
            if ($spanning_read_want{"$scaffold|$core_read_name"}) {
                print $sam_entry->get_original_line() . "\n";
            }
        }
    }
    
    {
        
        print STDERR "-outputting the spanning read info.\n";

        ## output the spanning read info
        my $spanning_read_info_file = "$bam_file.fusion_spanning_info";
        open (my $ofh, ">$spanning_read_info_file") or die "Error, cannot write to $spanning_read_info_file";
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

            print $ofh join("\t", $fusion_info, $num_spanning) . "\t" . join(",", @spanning_reads)
                . "\t$num_left_contrary_support\t$contrary_left_support_txt"
                . "\t$num_right_contrary_support\t$contrary_right_support_txt"
                . "\n";
        }
        close $ofh;
    }
    
    
    exit(0);
}

####
sub overlaps_exon {
    my ($genome_coords_aref, $exon_bounds_href) = @_;

    foreach my $coordset (@$genome_coords_aref) {
        my ($lend, $rend) = @$coordset;
        foreach my $exon_coordset (keys %$exon_bounds_href) {
            my ($e_lend, $e_rend) = split(/-/, $exon_coordset);
            if ($e_lend < $rend && $e_rend > $lend) {
                return(1);
            }
        }
    }
    return(0); # no dice
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
        $info =~ /gene_name \"([^\"]+)\"/ or die "Error, cannot parse gene_name from $info";
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


