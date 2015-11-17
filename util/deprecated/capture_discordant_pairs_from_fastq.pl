#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use SAM_reader;
use SAM_entry;
use Fastq_reader;


my $MIN_MAPPING_QUALITY = 10;


my $usage = "usage: $0 left.fq right.fq aligned.bam out_prefix\n\n";


my $left_fq = $ARGV[0] or die $usage;
my $right_fq = $ARGV[1] or die $usage;
my $aligned_bam = $ARGV[2] or die $usage;
my $out_prefix = $ARGV[3] or die $usage;


main: {

    my %fastq_reads_capture;

    my $c = 0;
    my $sam_reader = new SAM_reader($aligned_bam);
    while (my $sam_entry = $sam_reader->get_next()) {

        $c++;
        if ($c % 10000 == 0) {
            print STDERR "\r[$c]   ";
        }
        if ( (! $sam_entry->is_query_unmapped() )
             && 
             $sam_entry->get_mapping_quality() >= $MIN_MAPPING_QUALITY
             &&
             $sam_entry->is_mate_unmapped()
            ) {
            my $read_name = $sam_entry->get_read_name();
            $fastq_reads_capture{$read_name}++;
        }
    }
    
    my $num_reads_capture = scalar(keys %fastq_reads_capture);
    print STDERR "$0 : $num_reads_capture read pairs to capture as discordant.\n";
    
    my $new_left_fq = "${out_prefix}_1.fastq";
    my $new_right_fq = "${out_prefix}_2.fastq";
    
    
    foreach my $file_pair ( [$left_fq, $new_left_fq],
                            [$right_fq, $new_right_fq] ) {

        
        my ($in_file, $out_file) = @$file_pair;

        my $count_captured_reads = 0;

        open (my $ofh, ">$out_file") or die "Error, cannot write to $out_file";

        my $fq_reader = new Fastq_reader($in_file);
        
        while (my $fq_obj = $fq_reader->next()) {

            my $core_read_name = $fq_obj->get_core_read_name();
            
            if (exists $fastq_reads_capture{$core_read_name}) {
                
                print $ofh $fq_obj->get_fastq_record();
            
                $count_captured_reads++;

            }
        }
        print STDERR "$0 : $in_file : $count_captured_reads num_captured_reads\n";
        
        close $ofh;
        

        my $cmd = "gzip $out_file";
        my $ret = system($cmd);
        if ($ret) {
            die "Error, cmd: $cmd died with ret $ret";
        }
        
    }

    exit(0);
}
