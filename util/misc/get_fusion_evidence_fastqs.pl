#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");
use Fastq_reader;
use Process_cmd;
use DelimParser;

my $usage = "\n\n\tusage: $0 finspector.fusion_predictions.final  left.fq (right.fq|NONE) output_dir\n\n";

my $fusion_results_file = $ARGV[0] or die $usage;
my $left_fq = $ARGV[1] or die $usage;
my $right_fq = $ARGV[2] or die $usage;
my $outdir_name = $ARGV[3] or die $usage;


main: {

    ## get the core fragment names:
    
    my %core_frag_name_to_fusion_name;
    


    open (my $fh, $fusion_results_file) or die "Error, cannot open file $fusion_results_file";

    my $tab_reader = new DelimParser::Reader($fh, "\t");

    while (my $row = $tab_reader->get_row()) {
        my $fusion_name = $row->{'#FusionName'};
        my $junction_reads_list_txt = $row->{JunctionReads};
        my $spanning_frags_list_txt = $row->{SpanningFrags};
        
        my %reads = &parse_core_frag_names(join(",", $junction_reads_list_txt, $spanning_frags_list_txt));
        
        &append_reads_to_fusion($fusion_name, \%core_frag_name_to_fusion_name, \%reads);

    }
    close $fh;

    unless (-d $outdir_name) {
        &process_cmd("mkdir -p $outdir_name");
    }
    
    &write_fastq_files("_1", $left_fq, $outdir_name, \%core_frag_name_to_fusion_name);
    
    &write_fastq_files("_2", $right_fq, $outdir_name, \%core_frag_name_to_fusion_name) unless ($right_fq eq "NONE");
    

    print STDERR "\nDone.\n\n";
    
    exit(0);
    
}


####
sub write_fastq_files {
    my ($file_suffix, $fastq_file, $outdir_name, $core_frag_name_to_fusion_name_href) = @_;
    
    my %fusion_to_ofh;

    my $fastq_reader = new Fastq_reader($fastq_file);
    
    while (my $fq_record = $fastq_reader->next()) {
        
        my $core_read_name = $fq_record->get_core_read_name();
        #print "[$core_read_name]\n";
        
        if (my $fusion_name = $core_frag_name_to_fusion_name_href->{$core_read_name}) {

            my $record_text = $fq_record->get_fastq_record();
            
            my $ofh = $fusion_to_ofh{$fusion_name};
            unless ($ofh) {
                my $outfilename = join("/", $outdir_name, "${fusion_name}$file_suffix.fastq");
                print STDERR "-writing to $outfilename\n";
                open ($ofh, ">$outfilename") or die "Error, cannot write to file $outfilename";
                $fusion_to_ofh{$fusion_name} = $ofh;
            }
            print $ofh $record_text;
        }
    }
    foreach my $ofh (values %fusion_to_ofh) {
        close $ofh;
    }

    print STDERR "\nDone writing to $file_suffix fastq file series.\n\n";
    
    return;
}


####
sub append_reads_to_fusion {
    my ($fusion_name, $core_frag_name_to_fusion_name_href, $reads_href) = @_;

    foreach my $frag_name (keys %$reads_href) {
        $core_frag_name_to_fusion_name_href->{$frag_name} = $fusion_name ;
    }
    
    return;
}

####
sub parse_core_frag_names {
    my ($comma_delim_read_list_txt) = @_;

    my %core_frag_names;

    my @read_names = split(/,/, $comma_delim_read_list_txt);
    foreach my $read_name (@read_names) {
        $read_name =~ s|/[12]$||;
        $core_frag_names{$read_name} = 1;
    }

    
    return(%core_frag_names);
}


