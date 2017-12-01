#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use DelimParser;
use Data::Dumper;

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
#  --min_spanning_frags_only <int>        minimum number of rna-seq fragments required as fusion evidence if
#                                          there are no junction reads
#
#  --fusion_preds <string>               fusion predictions in star-F-ftm
#
#  Optional:
#
#  --require_LDAS  0|1>                 require long double anchor support for fusion-breakpoint reads
#                                        when there are no spanning frags.  default: 1 (true)
# 
#########################################################################################################



__EOUSAGE__

    ;



my $help_flag;

my $min_junction_reads;
my $min_sum_frags;
my $min_novel_junction_support;
my $min_spanning_frags_only;
my $fusion_preds;
my $REQUIRE_LDAS = 1;


&GetOptions ( 'h' => \$help_flag,
              'min_junction_reads=i' => \$min_junction_reads,
              'min_sum_frags=i' => \$min_sum_frags,
              'min_novel_junction_support=i' => \$min_novel_junction_support,
              'min_spanning_frags_only=i' => \$min_spanning_frags_only,
              'fusion_preds=s' => \$fusion_preds,
              'require_LDAS=i' => \$REQUIRE_LDAS,
    );



foreach my $var ($min_junction_reads, $min_sum_frags, 
                 $min_novel_junction_support, 
                 $min_spanning_frags_only,
                 $fusion_preds) {
    unless (defined $var) {
        die $usage;
    }
}


main: {
    
    
    open (my $fh, $fusion_preds) or die $!;
    my $tab_reader = new DelimParser::Reader($fh, "\t");
    my @column_headers = $tab_reader->get_column_headers();
    
    my $tab_writer = new DelimParser::Writer(*STDOUT, "\t", \@column_headers);
    
    while (my $row = $tab_reader->get_row()) {
        
        my $J = $row->{JunctionReadCount};
        my $S = $row->{SpanningFragCount};
                
        unless (defined($J)) {
            die "Error, JunctionReadCount not defined for entry: " . Dumper($row);
        }
        unless (defined($S)) {
            die "Error, SpanningFragCount not defined for entry: " . Dumper($row);
        }
                
        if ($J < $min_junction_reads) {
            next; 
        }
        
        my $sum_JS = $J + $S;
        if ($sum_JS < $min_sum_frags) {
            next;
        }
        
        if ($J == 0 && $S < $min_spanning_frags_only) {
            next;
        }
        
        my $splice_type = $row->{SpliceType};
        if ($splice_type eq 'INCL_NON_REF_SPLICE' && $J < $min_novel_junction_support) {
            next;
        }
        
        # require big anchors when no spanning support exists.
        my $large_breakpoint_anchored = $row->{LargeAnchorSupport};
        if ($S == 0 && $REQUIRE_LDAS && $large_breakpoint_anchored !~ /YES/) {
            next;
        }
        
        #######################################
        # if got here, passed all requirements.
        
        $tab_writer->write_row($row);
        
    }
    close $fh;
    
    
    exit(0);
}


        
