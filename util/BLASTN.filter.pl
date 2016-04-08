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
#  --genome_lib_dir <string>      genome lib dir (see FusionFilter.github.io)
#  --out_prefix <string>          output prefix for STAR-Fusion.filter  (adds .final and .final.abridged as outputs)
#  --max_promiscuity <int>               maximum number of partners allowed for a given fusion. Default: 3
#  -E <float>                     E-value threshold for blast searches (default: 0.001)
#
#
#
########################################################################  


__EOUSAGE__

    ;

my $help_flag;

my $fusion_preds_file;
my $out_prefix;
my $genome_lib_dir;

my $Evalue = 1e-3;
my $max_promiscuity = 3;


&GetOptions ( 'h' => \$help_flag, 
              
              'fusion_preds=s' => \$fusion_preds_file,
              'out_prefix=s' => \$out_prefix,
              'E=f' => \$Evalue,
              'max_promiscuity=i' => \$max_promiscuity,
              'genome_lib_dir=s' => \$genome_lib_dir,

    );


if ($help_flag) {
    die $usage;
}

unless ($fusion_preds_file && $out_prefix && defined($Evalue) && $max_promiscuity && $genome_lib_dir) {
    die $usage;
}


main: {

    my $star_fusion_fmt_file = "$fusion_preds_file.starFfmt";

    open (my $fh, $fusion_preds_file) or die "Error, cannot open file $fusion_preds_file";
    open (my $ofh, ">$star_fusion_fmt_file") or die  "Error, cannot write to $star_fusion_fmt_file";
    
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
    
    my $tab_writer = new DelimParser::Writer($ofh, "\t", \@adj_column_headers);
    
    while (my $row = $tab_reader->get_row()) {
        
        my $fusion_name = join("--", $row->{LeftGene}, $row->{RightGene});
        $row->{'#FusionName'} = $fusion_name;
        
        $tab_writer->write_row($row);
    }
    
    close $fh;
    close $ofh;
    
    ## now do the homology filter
    
    my $cmd = "$FindBin::Bin/../FusionFilter/blast_and_promiscuity_filter.pl " 
        . " --fusion_preds $star_fusion_fmt_file  "
        . " --max_promiscuity $max_promiscuity "
        . " --out_prefix $out_prefix "
        . " -E $Evalue "
        . " --genome_lib_dir $genome_lib_dir ";
    
    

    &process_cmd($cmd);

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
    
        
