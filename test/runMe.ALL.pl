#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Pipeliner;
use __GLOBALS__;

my $left_fq = "reads.left.simPE.fq.gz";
my $right_fq = "reads.right.simPE.fq.gz";

my $INSTALL_DIR = "$FindBin::Bin/../";

main: {

    
    my $pipeliner = new Pipeliner(-verbose => 1);
    
    ####################
    ## FusionInspector #
    ####################
    
    my $cmd = "$INSTALL_DIR/FusionInspector --fusions test_fusions.list,test_fusions.list2,test_fusions.list3 --gtf $FUSION_ANNOTATOR_LIB/gencode.v19.annotation.gtf.exons --genome_fa $FUSION_ANNOTATOR_LIB/Hg19.fa --cdna_fa $FUSION_ANNOTATOR_LIB/gencode.v19.annotation.gtf.exons.cdna --left_fq $left_fq --right $right_fq --out_dir Fusion_Inspector_ALL/ --out_prefix finspector --include_whole_genome --align_utils STAR,HISAT,GSNAP  ";
    
    if (@ARGV) {
        $cmd .= " --include_Trinity"
    }

    $pipeliner->add_commands( new Command($cmd, "Fusion_Inspector.ok") );
    
    ## Execute pipeline
    
    $pipeliner->run();

    
    
    exit(0);

}


                             



