#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Pipeliner;
use __GLOBALS__;

my $usage = "usage: $0 align_method=(STAR|HISAT|GSNAP) [INCLUDE_TRINITY_FLAG]\n\n";

my $method = $ARGV[0] or die $usage;
my $INCLUDE_TRINITY_FLAG = $ARGV[1] || 0;

unless ($method =~ /GSNAP|HISAT|STAR/) {
    die $usage;
}


my $left_fq = "reads.left.simPE.fq.gz";
my $right_fq = "reads.right.simPE.fq.gz";

my $INSTALL_DIR = "$FindBin::Bin/../";

unless ($ENV{CTAT_GENOME_LIB}) {
    die "Error, must set env var for CTAT_GENOME_LIB";
}

main: {

    
    my $pipeliner = new Pipeliner(-verbose => 1);
    
    ####################
    ## FusionInspector #
    ####################
    
    my $cmd = "$INSTALL_DIR/FusionInspector --fusions test_fusions.list,test_fusions.list2,test_fusions.list3 --genome_lib $ENV{CTAT_GENOME_LIB} --left_fq $left_fq --right $right_fq --out_dir Fusion_Inspector/ --out_prefix finspector --include_whole_genome --align_utils $method ";
    
    if ($INCLUDE_TRINITY_FLAG) {
        $cmd .= " --include_Trinity"
    }

    $pipeliner->add_commands( new Command($cmd, "Fusion_Inspector.ok") );
    
    ## Execute pipeline
    
    $pipeliner->run();

    
    
    exit(0);

}


                             



