#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Process_cmd;
use __GLOBALS__;


my $fusion_files_list = "fusion_targets.A.txt,fusion_targets.B.txt,fusion_targets.C.txt";

my $INSTALL_DIR = "$FindBin::Bin/../";

unless ($ENV{CTAT_GENOME_LIB}) {
    die "Error, must set env var for CTAT_GENOME_LIB";
}

main: {

    ####################
    ## FusionInspector #
    ####################
    
    my $cmd = "$INSTALL_DIR/FusionInspector --fusions $fusion_files_list --genome_lib $ENV{CTAT_GENOME_LIB}  --out_prefix finspector --prep_for_IGV  @ARGV";
    
    &process_cmd($cmd);
    
    exit(0);

}


                             



