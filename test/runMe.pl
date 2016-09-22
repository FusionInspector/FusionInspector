#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Pipeliner;
use __GLOBALS__;

my $usage = "usage: $0 align_method=STAR [INCLUDE_TRINITY_FLAG]\n\n";

my $method = $ARGV[0] || "STAR";
my $INCLUDE_TRINITY_FLAG = $ARGV[1] || 0;

unless ($method =~ /GSNAP|HISAT|STAR|orig/) {
    die $usage;
}


my $left_fq = "test.reads_1.fastq.gz";
my $right_fq = "test.reads_2.fastq.gz";

my $fusion_files_list = "fusion_targets.A.txt,fusion_targets.B.txt,fusion_targets.C.txt";

my $INSTALL_DIR = "$FindBin::Bin/../";

unless ($ENV{CTAT_GENOME_LIB}) {
    die "Error, must set env var for CTAT_GENOME_LIB";
}

main: {

    
    my $pipeliner = new Pipeliner(-verbose => 2);
    
    ####################
    ## FusionInspector #
    ####################

    my $outdir = "Fusion_Inspector-" . join("-", split(/,/, $method));
    
    my $cmd = "$INSTALL_DIR/FusionInspector --fusions $fusion_files_list --genome_lib $ENV{CTAT_GENOME_LIB} --left_fq $left_fq --right $right_fq --out_dir $outdir --out_prefix finspector --align_utils $method --prep_for_IGV --no_cleanup ";
    
    if ($INCLUDE_TRINITY_FLAG) {
        $cmd .= " --include_Trinity"
    }

    $pipeliner->add_commands( new Command($cmd, "$outdir.ok") );
    
    ## Execute pipeline
    
    $pipeliner->run();

    
    
    exit(0);

}


                             



