#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "\n\n\tusage: $0 left.fq finspector.fusion_predictions.final.abridged\n\n";

my $fq_filename = $ARGV[0] or die $usage;
my $finspector_results = $ARGV[1] or die $usage;

main: {

    my $num_frags = &get_num_total_frags($fq_filename);
    
    open (my $fh, $finspector_results) or die "Error, cannot open file $finspector_results";
    my $header = <$fh>;
    my @header_pts = split(/\t/, $header);
    
    $header = join("\t", @header_pts[0..2], "J_FFPM", "S_FFPM", @header_pts[3..$#header_pts]);
    print $header;

    while (<$fh>) {
        my @x = split(/\t/);
        my $J = $x[1];
        my $S = $x[2];

        my $J_FFPM = &compute_FFPM($J, $num_frags);
        my $S_FFPM = &compute_FFPM($S, $num_frags);

        print join("\t", @x[0..2], $J_FFPM, $S_FFPM, @x[3..$#x]);
    }
    close $fh;

    exit(0);

}

####
sub get_num_total_frags {
    my ($fq_file) = @_;

    my $num_lines = `cat $fq_file | wc -l`;
    $num_lines =~ /^(\d+)/ or die "Error, cannot extract line count from [$num_lines]";
    $num_lines = $1;

    my $num_seq_records = $num_lines / 4;

    return($num_seq_records);
}
            

####
sub compute_FFPM {
    my ($count_frags, $total_frags) = @_;

    my $ffpm = $count_frags / $total_frags * 1e6;

    $ffpm = sprintf("%.4f", $ffpm);
    
    return($ffpm);
}

