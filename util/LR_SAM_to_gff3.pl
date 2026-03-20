#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::RealBin/../PerlLib");
use SAM_reader;
use SAM_entry;
use Carp;
use List::Util qw(min max);
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);

my $sam_file;
my $help_flag;
my $allow_non_primary = 0;

&GetOptions(
    'help|h' => \$help_flag,
    'sam=s' => \$sam_file,
    'allow_non_primary' => \$allow_non_primary,
);

if ($help_flag || ! $sam_file) {
    die "usage: $0 --sam alignments.bam [--allow_non_primary]\n";
}

my %PATH_COUNTER;
my $sam_reader = new SAM_reader($sam_file);

while ($sam_reader->has_next()) {
    my $sam_entry = $sam_reader->get_next();

    if ((! is_primary_alignment($sam_entry)) && (! $allow_non_primary)) {
        next;
    }
    next if $sam_entry->is_query_unmapped();

    my $read_name = $sam_entry->get_read_name();
    next if $read_name =~ /\.p\d$/;

    my $sam_line = $sam_entry->get_original_line();
    my $NM = 0;
    if ($sam_line =~ /NM:i:(\d+)/) {
        $NM = $1;
    }
    else {
        die "Error, couldn't extract num mismatches from sam line: $sam_line";
    }

    my $cigar_align = $sam_entry->get_cigar_alignment();
    my $num_indel_nts = 0;
    while ($cigar_align =~ /(\d+)[DI]/g) {
        $num_indel_nts += $1;
    }
    my $num_mismatches = $NM - $num_indel_nts;
    confess "negative mismatch count" if $num_mismatches < 0;

    my $scaff_name = $sam_entry->get_scaffold_name();
    my $strand = $sam_entry->get_query_strand();
    my ($genome_coords_aref, $query_coords_aref) = $sam_entry->get_alignment_coords();

    my $align_len = 0;
    foreach my $coordset (@$genome_coords_aref) {
        my $seglen = abs($coordset->[1] - $coordset->[0]) + 1;
        $align_len += $seglen;
    }
    my $per_id = sprintf("%.1f", 100 - $num_mismatches / $align_len * 100);
    my $align_counter = "$read_name.p" . ++$PATH_COUNTER{$read_name};

    my @genome_n_trans_coords;
    while (@$genome_coords_aref) {
        my $genome_coordset_aref = shift @$genome_coords_aref;
        my $trans_coordset_aref = shift @$query_coords_aref;
        my ($genome_lend, $genome_rend) = @$genome_coordset_aref;
        my ($trans_lend, $trans_rend) = sort {$a <=> $b} @$trans_coordset_aref;
        push @genome_n_trans_coords, [$genome_lend, $genome_rend, $trans_lend, $trans_rend];
    }

    my @merged_coords;
    push @merged_coords, shift @genome_n_trans_coords;
    my $MERGE_DIST = 10;
    while (@genome_n_trans_coords) {
        my $coordset_ref = shift @genome_n_trans_coords;
        my $last_coordset_ref = $merged_coords[$#merged_coords];
        if ($coordset_ref->[0] - $last_coordset_ref->[1] <= $MERGE_DIST) {
            $last_coordset_ref->[1] = $coordset_ref->[1];
            if ($strand eq "+") {
                $last_coordset_ref->[3] = $coordset_ref->[3];
            }
            else {
                $last_coordset_ref->[2] = $coordset_ref->[2];
            }
        }
        else {
            push @merged_coords, $coordset_ref;
        }
    }

    foreach my $coordset_ref (@merged_coords) {
        my ($genome_lend, $genome_rend, $trans_lend, $trans_rend) = @$coordset_ref;
        print join("\t",
                   $scaff_name, "LR_SAM_to_gff3.pl", "match",
                   $genome_lend, $genome_rend, $per_id, $strand, ".",
                   "ID=$align_counter;Target=$read_name $trans_lend $trans_rend") . "\n";
    }
}

exit(0);

sub is_primary_alignment {
    my ($sam_entry) = @_;
    my $flag = $sam_entry->get_flag();

    return (($flag & 0x0100) == 0) ? 1 : 0;
}
