#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use DelimParser;
use FusionEM;

my $usage = "\n\tusage: $0 fusions.tsv\n\n";


my $fusions_tsv = $ARGV[0] or die $usage;


main: {
    
    open(my $fh, $fusions_tsv) or die "Error, cannot open file: $fusions_tsv";

    my $tab_reader = new DelimParser::Reader($fh, "\t");
    
    my $fusion_em = new FusionEM();
    
    my @fusion_rows;
    
    my @column_headers = $tab_reader->get_column_headers();

    while (my $row = $tab_reader->get_row()) {
        
        push (@fusion_rows, $row);
        
        my $fusion_name = $tab_reader->get_row_val($row, "#FusionName");
        
        my $left_breakpoint = $tab_reader->get_row_val($row, "LeftBreakpoint");
        my $right_breakpoint = $tab_reader->get_row_val($row, "RightBreakpoint");
        
        $fusion_name .= "::" . $left_breakpoint . "::" . $right_breakpoint;
        
        $row->{fusion_isoform_name} = $fusion_name; # save it for later...

        my $junction_reads_string = $tab_reader->get_row_val($row, "JunctionReads");
        my $spanning_frags_string = $tab_reader->get_row_val($row, "SpanningFrags");
       
        my @junction_reads;
        if ($junction_reads_string ne '.') {
            foreach my $junction_read (split(/,/, $junction_reads_string)) {
                $junction_read =~ s|/[12]$||;
                push (@junction_reads, $junction_read);
            }
        }

        my @spanning_frags;
        if ($spanning_frags_string ne '.') {
            push(@spanning_frags, split(/,/, $spanning_frags_string));
        }
        
        $fusion_em->add_fusion_transcript($fusion_name, \@junction_reads, \@spanning_frags);
        
    }

    $fusion_em->run();

    ## output with estimated J and S vals

    my @adjusted_column_headers = @column_headers;
    unless (grep { $_ eq "est_J" } @adjusted_column_headers) {
        @adjusted_column_headers = (@column_headers[0..2], "est_J", "est_S", @column_headers[3..$#column_headers]);
    }
    
    my $tab_writer = new DelimParser::Writer(*STDOUT, "\t", \@adjusted_column_headers);
    
    foreach my $row (@fusion_rows) {
        my $fusion_isoform_name = $row->{fusion_isoform_name};
        my ($est_J, $est_S) = $fusion_em->get_fusion_estimated_J_S($fusion_isoform_name);
        $row->{est_J} = sprintf("%.2f", $est_J);
        $row->{est_S} = sprintf("%.2f", $est_S);
        
        $tab_writer->write_row($row);
    }
    

    exit(0);
    
}


