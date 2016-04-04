#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Fasta_reader;
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


=incoming

0       #geneA
1       local_brkpt_A
2       chr_brkpt_A
3       geneB
4       local_brkpt_B
5       chr_brkpt_B
6       splice_type
7       junction_count
8       spanning_count
9       has_large_anchor_junction_support
10      junction_reads
11      spanning_reads
12      num_left_contrary_reads
13      left_contrary_reads
14      num_right_contrary_reads
15      right_contrary_reads
16      TAF_left
17      TAF_right
18      fusion_annotations

but want:

0       #fusion_name
1       JunctionReads
2       SpanningFrags
3       Splice_type
4       LeftGene
5       LeftBreakpoint
6       RightGene
7       RightBreakpoint
8       JunctionReads
9       SpanningFrags
10      annotations\tTrinityGG\t....

=cut


main: {

    my $star_fusion_fmt_file = "$fusion_preds_file.starFfmt";
    open (my $ofh, ">$star_fusion_fmt_file") or die  "Error, cannot write to $star_fusion_fmt_file";
    
    
    my @fusions;
    open (my $fh, $fusion_preds_file) or die "Error, cannot open file $fusion_preds_file";
    my $header = <$fh>;
    
    print $ofh join("\t", 
                    "#fusion_name", 
                    "JunctionReads", 
                    "SpanningFrags", 
                    "Splice_type", 
                    "LeftGene", 
                    "LeftBreakpoint", 
                    "RightGene", 
                    "RightBreakpoint", 
                    "JunctionReads", 
                    "SpanningFrags", 
                    "TAF_left",
                    "TAF_right",
                    "Annotations", 
                    "TrinityGG") . "\n";

    while (<$fh>) {
        if (/^\#/) { 
            next; 
        }
        chomp;
        my $line = $_;

        my ($geneA, 
            $local_chr_brkpt_A, 
            $chr_brkpt_A, 
            $geneB, 
            $local_chr_brkpt_B, 
            $chr_brkpt_B, 
            $splice_type, 
            $junction_count, 
            $spanning_count, 
            $has_lrg_anchor_support, 
            $junction_reads, 
            $spanning_reads, 
            $num_left_contrary_reads, 
            $left_contrary_reads, 
            $num_right_contrary_reads, 
            $right_contrary_reads, 
            $TAF_left, 
            $TAF_right, 
            @rest, # annotations, TrinityGG, ...
            ) = split(/\t/);
        
        my $fusion_name = "$geneA--$geneB";
        
        print $ofh join("\t", 
                        $fusion_name, 
                        $junction_count, 
                        $spanning_count, 
                        $splice_type, 
                        $geneA, 
                        $chr_brkpt_A, 
                        $geneB, 
                        $chr_brkpt_B, 
                        $junction_reads,
                        $spanning_reads, 
                        $TAF_left,
                        $TAF_right,
                        @rest) . "\n";
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
    
        
