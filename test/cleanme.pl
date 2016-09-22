#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;


## we delete all files we don't need in this directory. Be careful in case users try running it somewhere else, outside this dir.
chdir $FindBin::Bin or die "error, cannot cd to $FindBin::Bin";



my @files_to_keep = qw (cleanme.pl 
                        runMe.pl
runMe.hisat.pl
runMe.ALL.pl
README.txt

test.reads_1.fastq.gz
test.reads_2.fastq.gz
fusion_targets.A.txt
fusion_targets.B.txt
fusion_targets.C.txt

                        );


my %keep = map { + $_ => 1 } @files_to_keep;


`rm -rf ./Gsnap_Fusion` if (-d "Gsnap_Fusion");
`rm -rf ./Star_Fusion` if (-d "Star_Fusion");
`rm -rf ./Trinity_Fusion` if (-d "Trinity_Fusion");
`rm -rf ./Fusion_Inspector` if (-d "Fusion_Inspector");
`rm -rf ./trinity_out_dir` if (-d "trinity_out_dir");
`rm -rf ./prada_outdir` if (-d "prada_outdir");
`rm -rf ./soapfuse_outdir` if (-d "soapfuse_outdir");
`rm -rf ./Star_FGene` if (-d "Star_FGene");
`rm -rf ./_STAR*`;
`rm -rf ./Fusion_Inspector_ALL`;
`rm -rf ./Fusion_Inspector-*`;

foreach my $file (<*>) {
	
	if (-f $file && ! $keep{$file}) {
		print STDERR "-removing file: $file\n";
		unlink($file);
	}
}


exit(0);
