#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use DelimParser;


my $usage = "usage: $0 fusions.summary reads.frag_coords\n\n";

my $fusions_summary_file = $ARGV[0] or die $usage;
my $reads_frag_file = $ARGV[1] or die $usage;

## fusion summary: fi_test.fusion_preds.coalesced.summary

## frag coords: fi_test.fusion_preds.coalesced.summary.consolidated.cSorted.bam.frag_coords

main: {

    my %fusion_info;

    ## get junction reads counts and identification of spanning frags
    
    my %spanning_want;

    {
        
        open (my $fh, $fusions_summary_file) or die "Error, cannot open file $fusions_summary_file";
        
        my $tab_reader = new DelimParser::Reader($fh, "\t");
        
        while (my $row = $tab_reader->get_row()) {
                        
            my $geneA = $row->{LeftGene};
            my $geneB = $row->{RightGene};
            
            my $break_left = $row->{LeftLocalBreakpoint};
            my $break_right = $row->{RightLocalBreakpoint};

            my $num_junction_reads = $row->{JunctionReadCount};

            my $num_spanning_frags = $row->{SpanningFragCount};
            my @spanning_frags = split(/,/, $row->{SpanningFrags});
            
            my $fusion_name = $row->{'#FusionName'};
            
            foreach my $spanning_frag (@spanning_frags) {
                
                my $spanning_read_token = join("$;", $fusion_name, $spanning_frag);
                
                $spanning_want{$spanning_read_token} = undef; # just adding the key, set coords later.
            }
            
            push (@{$fusion_info{$fusion_name}}, { break_left => $break_left,
                                                   break_right => $break_right,

                                                   num_junction_reads => $num_junction_reads,
                                                   num_spanning_frags => $num_spanning_frags,
                                                   
                                                   spanning_frags => \@spanning_frags,
                                                   
                  } );
            
        }
    }

    ## get spanning_read_coords
    {
        open (my $fh, $reads_frag_file) or die "Error, cannot open file $reads_frag_file";
        while (<$fh>) {
            chomp;
            my ($fusion_contig, $frag_name, $lend, $rend) = split(/\t/);
            
            my $frag_token = join("$;", $fusion_contig, $frag_name);
            if (exists $spanning_want{$frag_token}) {
                $spanning_want{$frag_token} = "$lend-$rend";
            }
        }
        close $fh;
    }

    ## output JuncSpan fmt file
    
    # print header
    print join("\t", "#scaffold", "fusion_break_name", "break_left", "break_right",
               "num_junction_reads", "num_spanning_frags", "spanning_frag_coords") . "\n";

    foreach my $fusion_contig (sort keys %fusion_info) {
        
        my @breakpoint_structs = @{$fusion_info{$fusion_contig}};
        
        @breakpoint_structs = sort {$a->{break_left} <=> $b->{break_left}} @breakpoint_structs;

        foreach my $struct (@breakpoint_structs) {
            
            my $break_left = $struct->{break_left};
            my $break_right = $struct->{break_right};
            my $num_junction_reads = $struct->{num_junction_reads};
            my $num_spanning_frags = $struct->{num_spanning_frags};
            
            my @spanning_frags = @{$struct->{spanning_frags}};

            my $fusion_break_name = "$fusion_contig|$break_left-$break_right";

            my $outline = join("\t", $fusion_contig, $fusion_break_name, $break_left, $break_right, 
                               $num_junction_reads, $num_spanning_frags);

            my @span_coordsets;
            foreach my $spanning_frag (@spanning_frags) {
                my $frag_token = join("$;", $fusion_contig, $spanning_frag);
                if (my $coords = $spanning_want{$frag_token}) {
                    push (@span_coordsets, $coords);
                }
                else {
                    print STDERR "ERROR - no coordinates for $fusion_contig frag [$spanning_frag]\n";
                }
            }
            
            $outline .= "\t" . join(",", @span_coordsets);

            print "$outline\n";
        }
    }
    
    exit(0);
}


            
    

    
