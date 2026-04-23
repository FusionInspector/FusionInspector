#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

my $left_fq;
my $right_fq;
my $fusion;
my $support_type = "S";
my $out_prefix;
my $read_names_file;

my $usage = <<__EOUSAGE__;
usage: $0 --left left.fq[.gz] --right right.fq[.gz] --fusion FUSION --support_type S|J --out_prefix out_prefix [--read_names_file file]

Extract paired FASTQ records whose FusionInspector test-data comment line is:
  +FUSION|SUPPORT_TYPE|READ_CORE_NAME

If --read_names_file is provided, extract only those core read names instead.
__EOUSAGE__

GetOptions(
    "left=s" => \$left_fq,
    "right=s" => \$right_fq,
    "fusion=s" => \$fusion,
    "support_type=s" => \$support_type,
    "out_prefix=s" => \$out_prefix,
    "read_names_file=s" => \$read_names_file,
) or die $usage;

foreach my $required ($left_fq, $right_fq, $fusion, $support_type, $out_prefix) {
    die $usage unless defined $required && length $required;
}

main: {

    my $left_out = "${out_prefix}_1.fq";
    my $right_out = "${out_prefix}_2.fq";

    my %wanted_core_read = $read_names_file ? &load_read_names($read_names_file) : ();
    my $left_count = &extract_left_reads($left_fq, $left_out, \%wanted_core_read);
    my $right_count = &extract_right_reads($right_fq, $right_out, \%wanted_core_read);

    if ($left_count == 0) {
        die "Error, no left reads found for $fusion support type $support_type\n";
    }
    if ($left_count != $right_count) {
        die "Error, extracted $left_count left reads but $right_count right reads for $fusion support type $support_type\n";
    }

    print STDERR "Extracted $left_count read pairs for $fusion support type $support_type\n";

    exit(0);
}

sub extract_left_reads {
    my ($fq_file, $out_file, $wanted_core_read_href) = @_;

    my $fq_fh = &open_fastq($fq_file);
    open(my $out_fh, ">", $out_file) or die "Error, cannot write to $out_file";

    my $count = 0;
    while (my $record = &read_fastq_record($fq_fh)) {
        my ($header, $seq, $comment, $qual) = @$record;

        my $core_read;
        if ($read_names_file) {
            $core_read = &parse_core_read_from_header($header);
            next unless $wanted_core_read_href->{$core_read};
        }
        else {
            $core_read = &parse_supported_read($comment);
            next unless defined $core_read;
            print $out_fh join("", $header, $seq, $comment, $qual);
            $wanted_core_read_href->{$core_read} = 1;
            $count++;
            next;
        }

        print $out_fh join("", $header, $seq, $comment, $qual);
        $count++;
    }

    close $out_fh;
    close $fq_fh;

    return $count;
}

sub load_read_names {
    my ($filename) = @_;

    open(my $fh, $filename) or die "Error, cannot read $filename";

    my %read_names;
    while (my $line = <$fh>) {
        chomp $line;
        foreach my $read_name (split(/,/, $line)) {
            $read_name =~ s/^\s+|\s+$//g;
            next unless length $read_name;
            $read_name =~ s/^\&[^\@]+\@//;
            $read_name =~ s/\/[12]$//;
            $read_names{$read_name} = 1;
        }
    }

    close $fh;

    return %read_names;
}

sub extract_right_reads {
    my ($fq_file, $out_file, $wanted_core_read_href) = @_;

    my $fq_fh = &open_fastq($fq_file);
    open(my $out_fh, ">", $out_file) or die "Error, cannot write to $out_file";

    my $count = 0;
    while (my $record = &read_fastq_record($fq_fh)) {
        my ($header, $seq, $comment, $qual) = @$record;
        my $core_read = &parse_core_read_from_header($header);
        if ($wanted_core_read_href->{$core_read}) {
            print $out_fh join("", $header, $seq, $comment, $qual);
            $count++;
        }
    }

    close $out_fh;
    close $fq_fh;

    return $count;
}

sub parse_supported_read {
    my ($comment) = @_;

    chomp $comment;
    my $prefix = "+${fusion}|${support_type}|";

    if (index($comment, $prefix) == 0) {
        my $core_read = substr($comment, length($prefix));
        return $core_read if length $core_read;
    }

    return;
}

sub parse_core_read_from_header {
    my ($header) = @_;

    chomp $header;
    $header =~ s/^\@//;
    $header =~ s/\/[12]$//;

    return $header;
}

sub open_fastq {
    my ($fq_file) = @_;

    my $fh;
    if ($fq_file =~ /\.gz$/) {
        open($fh, "gunzip -c $fq_file |") or die "Error, cannot read $fq_file";
    }
    else {
        open($fh, $fq_file) or die "Error, cannot read $fq_file";
    }

    return $fh;
}

sub read_fastq_record {
    my ($fh) = @_;

    my $header = <$fh>;
    return unless defined $header;

    my $seq = <$fh>;
    my $comment = <$fh>;
    my $qual = <$fh>;

    unless (defined $seq && defined $comment && defined $qual) {
        die "Error, malformed FASTQ record at $header";
    }

    return [$header, $seq, $comment, $qual];
}
