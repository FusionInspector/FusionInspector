#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Fasta_reader;
use DelimParser;
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);                                                 


my $usage = <<__EOUSAGE__;

########################################################################
#
# Required:
#
#  --fusion_preds <string>        preliminary fusion prediction (file: finspector.fusion_preds.coalesced.summary)s
#
########################################################################  


__EOUSAGE__

    ;

my $help_flag;

my $fusion_preds_file;


&GetOptions ( 'h' => \$help_flag, 
              
              'fusion_preds=s' => \$fusion_preds_file,
    );


if ($help_flag) {
    die $usage;
}

unless ($fusion_preds_file) { 
    die $usage;
}


main: {

    open (my $fh, $fusion_preds_file) or die "Error, cannot open file $fusion_preds_file";
    
    my $tab_reader = new DelimParser::Reader($fh, "\t");
        
    my @column_headers = $tab_reader->get_column_headers();
    
    # reorder
    my %remove = ( JunctionReadCount => 1,
                   SpanningFragCount => 1 );
    
    my @adj_column_headers;
    foreach my $header (@column_headers) {
        if (! exists $remove{$header}) {
            push (@adj_column_headers, $header);
        }
    }

    @adj_column_headers = ("#FusionName", "JunctionReadCount", "SpanningFragCount", @adj_column_headers); # retain ordering for the rest.
    
    my $tab_writer = new DelimParser::Writer(*STDOUT, "\t", \@adj_column_headers);
    
    while (my $row = $tab_reader->get_row()) {

        my $left_gene = $row->{LeftGene};
        $left_gene =~ s/\^.*$//;  # just get the gene symbol if it's a compound gene_sym^gene_id

        my $right_gene = $row->{RightGene};
        $right_gene =~ s/\^.*$//;
        
        
        my $fusion_name = join("--", $left_gene, $right_gene);
        $row->{'#FusionName'} = $fusion_name;
        
        $tab_writer->write_row($row);
    }
    
    close $fh;
    
    exit(0);
        
}



####
sub process_cmd {
    my ($cmd) = @_;

    print STDERR "CMD: $cmd\n";
        
    my $ret = system($cmd);
    if ($ret) {

        die "Error, cmd $cmd died with ret $ret";
    }
    
    return;
}
    
        
