#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);


my $usage = <<__EOUSAGE__;

#########################################################################################################
#
#  --min_junction_reads <int>            minimum number of junction-spanning reads required
#
#  --min_sum_frags <int>                 minimum fusion support = ( # junction_reads + # spanning_frags )
#
#  --min_novel_junction_support <int>    minimum number of junction reads required if breakpoint 
#                                          lacks involvement of only reference junctions
#
#  --fusion_preds <string>               fusion predictions in star-F-ftm
#
#########################################################################################################



__EOUSAGE__

    ;



my $help_flag;

my $min_junction_reads;
my $min_sum_frags;
my $min_novel_junction_support;
my $fusion_preds;

&GetOptions ( 'h' => \$help_flag,
              'min_junction_reads=i' => \$min_junction_reads,
              'min_sum_frags=i' => \$min_sum_frags,
              'min_novel_junction_support=i' => \$min_novel_junction_support,
              'fusion_preds=s' => \$fusion_preds,
    );



foreach my $var ($min_junction_reads, $min_sum_frags, $min_novel_junction_support, $fusion_preds) {
    unless (defined $var) {
        die $usage;
    }
}


main: {
    
    
    open (my $fh, $fusion_preds) or die $!;
    while (<$fh>) {
        my $line = $_;
        if (/^\#/) { 
            print $line;
            next;
        }

        my @x = split(/\t/, $line);
        my $splice_type = $x[6];
        my $J = $x[7];
        my $S = $x[8];
        my $large_breakpoint_anchored = $x[9];

        if ($splice_type ne 'ONLY_REF_SPLICE') {
            unless ($J >= $min_novel_junction_support && $large_breakpoint_anchored eq 'YES') {
                next;
            }
        }

        # require big anchors when no spanning support exists.
        if ($S == 0 && $large_breakpoint_anchored eq 'NO') {
            next;
        }
        
        
        my $sum_JS = $J + $S;
        
        if ($J >= $min_junction_reads && $sum_JS >= $min_sum_frags) {

            print $line;
        }
    }
    close $fh;

    
    exit(0);
}


        
