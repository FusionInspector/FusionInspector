#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Process_cmd;
use __GLOBALS__;


my $left_fq = "test.reads_1.fastq.gz";
my $right_fq = "test.reads_2.fastq.gz";

my $fusion_files_list = "fusion_targets.A.txt,fusion_targets.B.txt,fusion_targets.C.txt";

my $INSTALL_DIR = "$FindBin::Bin/../";

unless ($ENV{CTAT_GENOME_LIB}) {
    die "Error, must set env var for CTAT_GENOME_LIB";
}

main: {

    ####################
    ## FusionInspector #
    ####################
    
    my $cmd = "$INSTALL_DIR/FusionInspector --fusions $fusion_files_list --genome_lib $ENV{CTAT_GENOME_LIB} --left_fq $left_fq --right $right_fq --out_prefix finspector --prep_for_IGV --no_cleanup @ARGV";
    
    &process_cmd($cmd);
    
    exit(0);

}


                             



