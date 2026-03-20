#!/usr/bin/env perl

use strict;
use warnings;
use Carp;

use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Overlap_piler;
use List::Util qw(min max);
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);

my ($FI_gtf_filename, $LR_gff3_filename, $seq_similar_gff3_filename, $SNAP_dist, $min_trans_overlap_length);
my $PSEUDOCOUNT = 1;
my $DEBUG = 0;

&GetOptions(
    'FI_gtf=s' => \$FI_gtf_filename,
    'LR_gff3=s' => \$LR_gff3_filename,
    'seq_similar_gff3=s' => \$seq_similar_gff3_filename,
    'snap_dist=i' => \$SNAP_dist,
    'min_trans_overlap_length=i' => \$min_trans_overlap_length,
    'debug|d' => \$DEBUG,
);

unless ($FI_gtf_filename && $LR_gff3_filename && $seq_similar_gff3_filename && defined($SNAP_dist) && defined($min_trans_overlap_length)) {
    die "usage: $0 --FI_gtf finspector.gtf --LR_gff3 alignments.gff3 --seq_similar_gff3 seqsimilar.gff3 --snap_dist 3 --min_trans_overlap_length 100\n";
}

my $DONOR_TYPE = "DONOR";
my $ACCEPTOR_TYPE = "ACCEPTOR";
my $NA_TYPE = "NA";

my %orig_coord_info;
my %scaffold_to_gene_coordsets;
my %scaffold_to_gene_trans_to_coordsets;
parse_FI_gtf_filename($FI_gtf_filename, \%orig_coord_info, \%scaffold_to_gene_coordsets, \%scaffold_to_gene_trans_to_coordsets);
my %scaffold_to_orig_coords = organize_original_coordinate_info(\%orig_coord_info);
my %seqsimilar_regions = parse_seqsimilar_gff3($seq_similar_gff3_filename);
my %scaffold_to_LR_coords = parse_LR_alignment_gff3_file($LR_gff3_filename);
my %scaffold_gene_bounds;

my %LR_fusion_trans_ids;

foreach my $scaffold (sort keys %scaffold_to_gene_coordsets) {
    my ($left_gene, $right_gene) = split(/--/, $scaffold);

    my ($geneA_coords_href, $geneB_coords_href) = get_gene_coords($scaffold, $scaffold_to_gene_coordsets{$scaffold});
    my ($transA_all_coords_aref, $transB_all_coords_aref) = get_trans_coordsets($scaffold, $scaffold_to_gene_trans_to_coordsets{$scaffold});

    my $seqsimilar_regions_aref = $seqsimilar_regions{$scaffold};
    if ($seqsimilar_regions_aref) {
        $transA_all_coords_aref = exclude_seqsimilar_regions($transA_all_coords_aref, $seqsimilar_regions_aref);
        $transB_all_coords_aref = exclude_seqsimilar_regions($transB_all_coords_aref, $seqsimilar_regions_aref);
    }
    next unless (@$transA_all_coords_aref && @$transB_all_coords_aref);

    $transA_all_coords_aref = collapse_overlapping_trans_segments($transA_all_coords_aref);
    $transB_all_coords_aref = collapse_overlapping_trans_segments($transB_all_coords_aref);

    my $geneA_max = max(keys %$geneA_coords_href);
    my $geneB_min = min(keys %$geneB_coords_href);
    $scaffold_gene_bounds{$scaffold} = [$geneA_max, $geneB_min];
    next unless exists $scaffold_to_LR_coords{$scaffold};

    foreach my $LR_acc (sort keys %{$scaffold_to_LR_coords{$scaffold}}) {
        my @LR_coordsets = sort {$a->[0] <=> $b->[0]} @{$scaffold_to_LR_coords{$scaffold}->{$LR_acc}};
        next if scalar(@LR_coordsets) < 2;

        my $min_LR_coord = $LR_coordsets[0]->[0];
        my $max_LR_coord = $LR_coordsets[$#LR_coordsets]->[1];
        next unless ($min_LR_coord < $geneA_max && $max_LR_coord > $geneB_min);

        my @left_gene_align_coords = grep { $_->[0] < $geneA_max } @LR_coordsets;
        my @right_gene_align_coords = grep { $_->[1] > $geneB_min } @LR_coordsets;

        next unless has_exon_overlapping_segment(\@left_gene_align_coords, $transA_all_coords_aref);
        next unless has_exon_overlapping_segment(\@right_gene_align_coords, $transB_all_coords_aref);

        my $left_gene_overlapped_bases = sum_overlaps(\@left_gene_align_coords, $transA_all_coords_aref);
        next if $left_gene_overlapped_bases < $min_trans_overlap_length;
        my $right_gene_overlapped_bases = sum_overlaps(\@right_gene_align_coords, $transB_all_coords_aref);
        next if $right_gene_overlapped_bases < $min_trans_overlap_length;

        my ($break_left, $break_right) = get_breakpoint_coords(\@LR_coordsets, $geneA_max, $geneB_min);
        $LR_fusion_trans_ids{$LR_acc}->{"$scaffold:$break_left-$break_right"} = 1;
    }
}

report_LR_fusions(\%LR_fusion_trans_ids, \%orig_coord_info, \%scaffold_to_orig_coords, \%scaffold_to_LR_coords, \%scaffold_gene_bounds);
exit(0);

sub report_LR_fusions {
    my ($LR_ids_href, $orig_coord_info_href, $scaffold_to_orig_coords_href, $scaffold_to_LR_coords_href, $scaffold_gene_bounds_href) = @_;
    my %scaff_breakpoint_to_read_support;
    foreach my $LR_id (sort keys %$LR_ids_href) {
        foreach my $scaff_breakpoint (sort keys %{$LR_ids_href->{$LR_id}}) {
            push @{$scaff_breakpoint_to_read_support{$scaff_breakpoint}}, $LR_id;
        }
    }

    my @fusion_structs;
    foreach my $breakpoint (sort keys %scaff_breakpoint_to_read_support) {
        my @LR_reads = @{$scaff_breakpoint_to_read_support{$breakpoint}};
        my $num_reads = scalar(@LR_reads);
        my ($scaffold, $breakpoint_coords) = split(/:/, $breakpoint);
        my ($break_lend, $break_rend) = split(/-/, $breakpoint_coords);

        my ($adj_break_lend, $left_genome_breakpoint, $left_ref_splice_mapping) =
            infer_genome_breakpoint_from_local_coord($break_lend, $orig_coord_info_href->{$scaffold}, $scaffold_to_orig_coords_href->{$scaffold}, $DONOR_TYPE, $SNAP_dist);
        $break_lend = $adj_break_lend;

        my ($adj_break_rend, $right_genome_breakpoint, $right_ref_splice_mapping) =
            infer_genome_breakpoint_from_local_coord($break_rend, $orig_coord_info_href->{$scaffold}, $scaffold_to_orig_coords_href->{$scaffold}, $ACCEPTOR_TYPE, $SNAP_dist);
        $break_rend = $adj_break_rend;

        my $splice_type = ($left_ref_splice_mapping == 1 && $right_ref_splice_mapping == 1) ? "ONLY_REF_SPLICE" : "INCL_NON_REF_SPLICE";
        my ($left_contrary_reads_aref, $right_contrary_reads_aref) =
            get_counter_fusion_reads($scaffold, $break_lend, $break_rend, $scaffold_to_LR_coords_href->{$scaffold}, $scaffold_gene_bounds_href->{$scaffold});
        my ($left_gene, $right_gene) = split(/--/, $scaffold);
        push @fusion_structs, {
            fusion_name => $scaffold,
            LeftGene => $left_gene,
            RightGene => $right_gene,
            LeftLocalBreakpoint => $break_lend,
            RightLocalBreakpoint => $break_rend,
            LeftBreakpoint => $left_genome_breakpoint,
            RightBreakpoint => $right_genome_breakpoint,
            JunctionReadCount => $num_reads,
            JunctionReads => \@LR_reads,
            SpliceType => $splice_type,
            CounterFusionLeftReads => $left_contrary_reads_aref,
            CounterFusionRightReads => $right_contrary_reads_aref,
        };
    }

    @fusion_structs = merge_identical_breakpoints(@fusion_structs);
    @fusion_structs = reverse sort {
        $a->{JunctionReadCount} <=> $b->{JunctionReadCount}
        ||
        $b->{fusion_name} cmp $a->{fusion_name}
    } @fusion_structs;

    print join("\t",
               "#FusionName", "JunctionReadCount", "SpanningFragCount", "est_J", "est_S",
               "LeftGene", "LeftLocalBreakpoint", "LeftBreakpoint",
               "RightGene", "RightLocalBreakpoint", "RightBreakpoint",
               "SpliceType", "LargeAnchorSupport",
               "JunctionReads", "SpanningFrags",
               "NumCounterFusionLeft", "CounterFusionLeftReads",
               "NumCounterFusionRight", "CounterFusionRightReads",
               "FAR_left", "FAR_right") . "\n";

    foreach my $fusion (@fusion_structs) {
        my $j = $fusion->{JunctionReadCount};
        my @left_counter_reads = @{$fusion->{CounterFusionLeftReads}};
        my @right_counter_reads = @{$fusion->{CounterFusionRightReads}};
        my $num_left_counter_reads = scalar(@left_counter_reads);
        my $num_right_counter_reads = scalar(@right_counter_reads);
        @left_counter_reads = (".") unless @left_counter_reads;
        @right_counter_reads = (".") unless @right_counter_reads;
        my $far_left = sprintf("%.2f", ($j + $PSEUDOCOUNT) / ($num_left_counter_reads + $PSEUDOCOUNT));
        my $far_right = sprintf("%.2f", ($j + $PSEUDOCOUNT) / ($num_right_counter_reads + $PSEUDOCOUNT));
        print join("\t",
                   $fusion->{fusion_name},
                   $j,
                   0,
                   sprintf("%.2f", $j),
                   sprintf("%.2f", 0),
                   $fusion->{LeftGene},
                   $fusion->{LeftLocalBreakpoint},
                   $fusion->{LeftBreakpoint},
                   $fusion->{RightGene},
                   $fusion->{RightLocalBreakpoint},
                   $fusion->{RightBreakpoint},
                   $fusion->{SpliceType},
                   "YES",
                   join(",", @{$fusion->{JunctionReads}}),
                   ".",
                   $num_left_counter_reads,
                   join(",", @left_counter_reads),
                   $num_right_counter_reads,
                   join(",", @right_counter_reads),
                   $far_left,
                   $far_right) . "\n";
    }
}

sub get_counter_fusion_reads {
    my ($scaffold, $break_lend, $break_rend, $scaffold_read_coords_href, $gene_bounds_aref) = @_;

    my @left_contrary_reads;
    my @right_contrary_reads;
    return (\@left_contrary_reads, \@right_contrary_reads) unless ($scaffold_read_coords_href && $gene_bounds_aref);

    my ($geneA_max, $geneB_min) = @$gene_bounds_aref;

    foreach my $LR_acc (sort keys %$scaffold_read_coords_href) {
        my @LR_coordsets = sort {$a->[0] <=> $b->[0]} @{$scaffold_read_coords_href->{$LR_acc}};
        next unless @LR_coordsets;

        my $min_LR_coord = $LR_coordsets[0]->[0];
        my $max_LR_coord = $LR_coordsets[$#LR_coordsets]->[1];

        if ($min_LR_coord < $break_lend && $break_lend < $max_LR_coord && $max_LR_coord < $geneB_min) {
            push @left_contrary_reads, $LR_acc;
        }
        elsif ($min_LR_coord > $geneA_max && $min_LR_coord < $break_rend && $break_rend < $max_LR_coord) {
            push @right_contrary_reads, $LR_acc;
        }
    }

    return (\@left_contrary_reads, \@right_contrary_reads);
}

sub get_gene_coords {
    my ($scaffold, $genes_to_coords_href) = @_;
    my @gene_coords_hrefs;
    foreach my $gene (split(/--/, $scaffold)) {
        my @coords = @{$genes_to_coords_href->{$gene}};
        my %gene_coords;
        foreach my $coordpair (@coords) {
            my ($lend, $rend) = @$coordpair;
            $gene_coords{$lend} = 1;
            $gene_coords{$rend} = 1;
        }
        push @gene_coords_hrefs, \%gene_coords;
    }
    return @gene_coords_hrefs;
}

sub parse_FI_gtf_filename {
    my ($FI_gtf_filename, $orig_coord_info_href, $scaff_to_gene_to_coords_href, $scaff_to_gene_trans_to_coords_href) = @_;
    open(my $fh, $FI_gtf_filename) or die "Error, cannot open file $FI_gtf_filename";
    while (<$fh>) {
        chomp;
        my @x = split(/\t/);
        next unless $x[2] eq 'exon';
        my $scaffold_id = $x[0];
        my $info = $x[8];
        my $gene_id = "";
        if ($info =~ /gene_name \"([^\"]+)\"/) {
            $gene_id = $1;
        }
        elsif ($info =~ /FI_gene_label \"([^\"]+)\"/) {
            $gene_id = $1;
            my @g = split(/\^/, $gene_id);
            $gene_id = $g[0];
        }
        else {
            die "Error extracting gene_name from $info";
        }
        $info =~ /transcript_id \"([^\"]+)\"/ or die "Error extracting transcript_id from $info";
        my $transcript_id = $1;
        my ($lend, $rend) = ($x[3], $x[4]);
        push @{$scaff_to_gene_to_coords_href->{$scaffold_id}->{$gene_id}}, [$lend, $rend];
        push @{$scaff_to_gene_trans_to_coords_href->{$scaffold_id}->{$gene_id}->{$transcript_id}}, [$lend, $rend];
        my $strand = $x[6];
        confess "FI contigs should be + strand" unless $strand eq '+';
        $info =~ /orig_coord_info \"([^,]+),(\d+),(\d+),([+-])\"/ or die "Error parsing orig_coord_info";
        my ($orig_chr, $orig_lend, $orig_rend, $orig_orient) = ($1, $2, $3, $4);
        my ($orig_end5, $orig_end3) = ($orig_orient eq '+') ? ($orig_lend, $orig_rend) : ($orig_rend, $orig_lend);
        $orig_coord_info_href->{$scaffold_id}->{$lend} = { chrom => $orig_chr, coord => $orig_end5, orient => $orig_orient, contig_coord => $lend, splice_junc => $NA_TYPE };
        $orig_coord_info_href->{$scaffold_id}->{$rend} = { chrom => $orig_chr, coord => $orig_end3, orient => $orig_orient, contig_coord => $rend, splice_junc => $NA_TYPE };
    }
    close $fh;
    foreach my $scaffold_id (sort keys %$scaff_to_gene_trans_to_coords_href) {
        foreach my $fusion_gene (split(/--/, $scaffold_id)) {
            foreach my $transcript_id (sort keys %{$scaff_to_gene_trans_to_coords_href->{$scaffold_id}->{$fusion_gene}}) {
                my @coordsets = sort {$a->[0] <=> $b->[0]} @{$scaff_to_gene_trans_to_coords_href->{$scaffold_id}->{$fusion_gene}->{$transcript_id}};
                for (my $i = 0; $i < $#coordsets; $i++) {
                    my $left_junc = $coordsets[$i]->[1];
                    my $right_junc = $coordsets[$i+1]->[0];
                    $orig_coord_info_href->{$scaffold_id}->{$left_junc}->{splice_junc} = $DONOR_TYPE;
                    $orig_coord_info_href->{$scaffold_id}->{$right_junc}->{splice_junc} = $ACCEPTOR_TYPE;
                }
            }
        }
    }
}

sub parse_LR_alignment_gff3_file {
    my ($LR_gff3_filename) = @_;
    my %scaffold_to_read_coords;
    open(my $fh, $LR_gff3_filename) or die $!;
    while (<$fh>) {
        next if /^\#/;
        next unless /\w/;
        chomp;
        my @x = split(/\t/);
        my ($scaff, $lend, $rend, $info) = ($x[0], $x[3], $x[4], $x[8]);
        my $LR_id;
        if ($info =~ /ID=([^;]+)\.p\d+;/) {
            $LR_id = $1;
        } else {
            die "Error, cannot find LR ID from $info";
        }
        push @{$scaffold_to_read_coords{$scaff}->{$LR_id}}, [$lend, $rend];
    }
    close $fh;
    return %scaffold_to_read_coords;
}

sub get_breakpoint_coords {
    my ($LR_coordsets_aref, $geneA_max, $geneB_min) = @_;
    for (my $i = 0; $i < $#$LR_coordsets_aref; $i++) {
        my $left_end = $LR_coordsets_aref->[$i]->[1];
        my $right_end = $LR_coordsets_aref->[$i+1]->[0];
        if (is_closer($left_end, $geneA_max, $geneB_min) && is_closer($right_end, $geneB_min, $geneA_max)) {
            return($left_end, $right_end);
        }
    }
    confess("Error, not finding proper fusion breakpoint");
}

sub is_closer {
    my ($coord_A, $coord_B, $coord_C) = @_;
    return abs($coord_A - $coord_B) < abs($coord_A - $coord_C);
}

sub organize_original_coordinate_info {
    my ($orig_coord_info_href) = @_;
    my %scaffold_to_orig_coords;
    foreach my $scaffold (sort keys %$orig_coord_info_href) {
        my @coordinates = sort {$a <=> $b} keys %{$orig_coord_info_href->{$scaffold}};
        my @coord_structs = map { $orig_coord_info_href->{$scaffold}->{$_} } @coordinates;
        @coord_structs = sort {$a->{contig_coord} <=> $b->{contig_coord}} @coord_structs;
        $scaffold_to_orig_coords{$scaffold} = \@coord_structs;
    }
    return %scaffold_to_orig_coords;
}

sub infer_genome_breakpoint_from_local_coord {
    my ($break_coord, $scaffold_orig_coord_info_href, $scaffold_coordinate_mappings_aref, $SPLICE_SITE_TYPE, $snap_dist) = @_;
    my $struct = $scaffold_orig_coord_info_href->{$break_coord};
    if (defined($struct) && $struct->{splice_junc} eq $SPLICE_SITE_TYPE) {
        return($break_coord, join(":", $struct->{chrom}, $struct->{coord}, $struct->{orient}), 1);
    }
    my @closest_structs;
    foreach my $s (@$scaffold_coordinate_mappings_aref) {
        my $delta = $break_coord - $s->{contig_coord};
        push @closest_structs, { struct => $s, delta => $delta, abs_delta => abs($delta) };
    }
    @closest_structs = sort {$a->{abs_delta} <=> $b->{abs_delta}} @closest_structs;
    my $closest_struct = shift @closest_structs;
    if ($closest_struct->{struct}->{splice_junc} ne $SPLICE_SITE_TYPE) {
        foreach my $other_struct (@closest_structs) {
            last if $other_struct->{abs_delta} > $snap_dist;
            if ($other_struct->{struct}->{splice_junc} eq $SPLICE_SITE_TYPE) {
                $closest_struct = $other_struct;
                last;
            }
        }
    }
    my ($chrom, $orient) = ($closest_struct->{struct}->{chrom}, $closest_struct->{struct}->{orient});
    my $coord;
    if ($closest_struct->{abs_delta} <= $snap_dist && $closest_struct->{struct}->{splice_junc} eq $SPLICE_SITE_TYPE) {
        $coord = $closest_struct->{struct}->{coord};
        return($closest_struct->{struct}->{contig_coord}, "$chrom:$coord:$orient", 1);
    }
    if ($orient eq '+') {
        $coord = $closest_struct->{struct}->{coord} + $closest_struct->{delta};
    }
    else {
        $coord = $closest_struct->{struct}->{coord} - $closest_struct->{delta};
    }
    return($break_coord, "$chrom:$coord:$orient", 0);
}

sub merge_identical_breakpoints {
    my (@fusion_structs) = @_;
    my %fusion_token_to_consolidated_fusions;
    foreach my $fusion (@fusion_structs) {
        my $token = join("$;", $fusion->{fusion_name}, $fusion->{LeftLocalBreakpoint}, $fusion->{RightLocalBreakpoint}, $fusion->{LeftBreakpoint}, $fusion->{RightBreakpoint}, $fusion->{SpliceType});
        if (my $existing = $fusion_token_to_consolidated_fusions{$token}) {
            $existing->{JunctionReadCount} += $fusion->{JunctionReadCount};
            push @{$existing->{JunctionReads}}, @{$fusion->{JunctionReads}};
        }
        else {
            $fusion_token_to_consolidated_fusions{$token} = $fusion;
        }
    }
    return values %fusion_token_to_consolidated_fusions;
}

sub get_trans_coordsets {
    my ($scaffold, $scaffold_to_trans_coords_href) = @_;
    my ($left_gene, $right_gene) = split(/--/, $scaffold);
    my @left_exon_coords = get_trans_coordsets_single_gene($scaffold_to_trans_coords_href->{$left_gene});
    my @right_exon_coords = get_trans_coordsets_single_gene($scaffold_to_trans_coords_href->{$right_gene});
    return(\@left_exon_coords, \@right_exon_coords);
}

sub get_trans_coordsets_single_gene {
    my ($trans_coords_href) = @_;
    my @all_coord_pairs;
    foreach my $trans_id (sort keys %$trans_coords_href) {
        foreach my $coordpair_href (values %$trans_coords_href) {
            push @all_coord_pairs, @$coordpair_href;
        }
    }
    return @all_coord_pairs;
}

sub has_exon_overlapping_segment {
    my ($LR_align_coords_aref, $all_trans_coords_aref) = @_;
    foreach my $align_coordset_aref (@$LR_align_coords_aref) {
        my ($align_lend, $align_rend) = sort {$a <=> $b} @$align_coordset_aref;
        foreach my $trans_coordset (@$all_trans_coords_aref) {
            my ($trans_lend, $trans_rend) = sort {$a <=> $b} @$trans_coordset;
            if ($trans_lend < $align_rend && $trans_rend > $align_lend) {
                return 1;
            }
        }
    }
    return 0;
}

sub parse_seqsimilar_gff3 {
    my ($seqsimilar_gff3_file) = @_;
    my %seqsimilar_regions;
    open(my $fh, $seqsimilar_gff3_file) or die "Error, cannot open file: $seqsimilar_gff3_file";
    while(<$fh>) {
        chomp;
        my @x = split(/\t/);
        if ($x[3] =~ /\d+/ && $x[4] =~ /\d+/) {
            push(@{$seqsimilar_regions{$x[0]}}, [$x[3], $x[4]]);
        }
    }
    close $fh;
    return %seqsimilar_regions;
}

sub exclude_seqsimilar_regions {
    my ($trans_aref, $seqsimilar_regions_aref) = @_;
    my @surviving_trans_coords;
    foreach my $trans_coordset (@$trans_aref) {
        my ($lend, $rend) = @$trans_coordset;
        my $found_overlap = 0;
        foreach my $seqsimilar_region (@$seqsimilar_regions_aref) {
            my ($region_lend, $region_rend) = @$seqsimilar_region;
            if ($lend <= $region_rend && $rend >= $region_lend) {
                $found_overlap = 1;
                last;
            }
        }
        push @surviving_trans_coords, [$lend, $rend] unless $found_overlap;
    }
    return \@surviving_trans_coords;
}

sub collapse_overlapping_trans_segments {
    my ($trans_coords_aref) = @_;
    my @coords = @$trans_coords_aref;
    return $trans_coords_aref if scalar(@coords) < 2;
    my @collapsed_coords = Overlap_piler::simple_coordsets_collapser(@coords);
    return \@collapsed_coords;
}

sub sum_overlaps {
    my ($coordsets_A_aref, $coordsets_B_aref) = @_;
    my $sum = 0;

    foreach my $coordset_A_aref (@$coordsets_A_aref) {
        foreach my $coordset_B_aref (@$coordsets_B_aref) {
            $sum += overlap_length($coordset_A_aref, $coordset_B_aref);
        }
    }

    return $sum;
}

sub overlap_length {
    my ($coordsA_aref, $coordsB_aref) = @_;
    return 0 unless overlaps($coordsA_aref, $coordsB_aref);

    my ($lendA, $rendA) = sort {$a <=> $b} @$coordsA_aref;
    my ($lendB, $rendB) = sort {$a <=> $b} @$coordsB_aref;

    my $overlap_lend = max($lendA, $lendB);
    my $overlap_rend = min($rendA, $rendB);

    return $overlap_rend - $overlap_lend + 1;
}

sub overlaps {
    my ($coordsA_aref, $coordsB_aref) = @_;

    my ($lendA, $rendA) = sort {$a <=> $b} @$coordsA_aref;
    my ($lendB, $rendB) = sort {$a <=> $b} @$coordsB_aref;

    return ($lendA <= $rendB && $rendA >= $lendB) ? 1 : 0;
}
