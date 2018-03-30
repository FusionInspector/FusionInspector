#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/../../PerlLib";
use Process_cmd;

my $usage = "\n\n\tusage: $0 samples.txt target_fusions.file\n\n";

my $samples_file = $ARGV[0] or die $usage;
my $target_fusions_file = $ARGV[1] or die $usage;

$target_fusions_file = &ensure_full_path($target_fusions_file);

main: {

    open(my $fh, $samples_file) or die "error, cannot open file $samples_file";
    while(<$fh>) {
        chomp;
        my ($sample_name, $left_fq, $right_fq) = split(/\s+/);

        my $cmd = "$FindBin::Bin/../../FusionInspector --fusions $target_fusions_file "
            . " --left_fq $left_fq --right_fq $right_fq --out_dir FI-$sample_name ";
        
        print "$cmd\n";

    }
    close $fh;

    exit(0);
}


        

 
    
