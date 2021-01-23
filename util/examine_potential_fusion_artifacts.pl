#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);


use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use DelimParser;
use Fasta_reader;
use Process_cmd;

my $UTILDIR = "$FindBin::Bin";


my $usage = <<__EOUSAGE__;

##############################################################################################
#
# --fusion_contigs_fasta <string>    finspector.fa
#
# --fusion_contigs_gtf <string>      finspector.gtf
#
# --output_directory <string>        /path/to/output/directory
#                                    (note, each fusion will be given a separate subdirectory)
#
# --fusion_results_tsv <string>      fusioninspector predictions (abridged.tsv file is fine)
# 
# --CPU <int>                        number of CPUs for parallel processing
#
# # optional
#
# --fusions_restrict <string>        file containing lists of fusions to restrict analysis to
#
###############################################################################################


__EOUSAGE__

    ;

my $help_flag;
my $fusion_contigs_fasta;
my $fusion_contigs_gtf;
my $output_directory;
my $fusion_results_tsv;
my $CPU;
my $fusions_restrict_file;


&GetOptions ( 'help|h' => \$help_flag,
              'fusion_contigs_fasta=s' => \$fusion_contigs_fasta,
              'fusion_contigs_gtf=s' => \$fusion_contigs_gtf,
              'output_directory=s' => \$output_directory,
              'fusion_results_tsv=s' => \$fusion_results_tsv,
              'CPU=i' => \$CPU,
              'fusions_restrict=s' => \$fusions_restrict_file,
              
    );

if ($help_flag) {
    die $usage; 
}

unless ($fusion_contigs_fasta && $fusion_contigs_gtf && $output_directory && $fusion_results_tsv && defined($CPU)) {
    die $usage;
}


my %FUSIONS_RESTRICT;
if ($fusions_restrict_file) {
    open(my $fh, $fusions_restrict_file) or die "Error, cannot open file $fusions_restrict_file";
    while(<$fh>) {
        chomp;
        $FUSIONS_RESTRICT{$_} = 1;
    }
}


main: {

    
    if (! -d $output_directory) {
        process_cmd("mkdir -p $output_directory");
    }
        
    ## Get fusions predictions:
    
    print STDERR "-parsing fusions: $fusion_results_tsv\n";
    my %fusion_name_to_rows;
    my @column_headers;
    {
        open(my $fh, $fusion_results_tsv) or die "Error, cannot open file: $fusion_results_tsv";
        my $delim_reader = new DelimParser::Reader($fh, "\t");
        @column_headers = $delim_reader->get_column_headers();
        
        
        while (my $row = $delim_reader->get_row()) {
            my $fusion_name = (exists $row->{"#FusionName"}) ? $delim_reader->get_row_val($row, "#FusionName")
                : (exists $row->{fusion_name}) ? $delim_reader->get_row_val($row, "fusion_name") 
                : confess "Error, cannot find #FusionName or fusion_name in [@column_headers] column headers";
            
            push(@{$fusion_name_to_rows{$fusion_name}}, $row);
        }
    }
    
    ## Get genome fasta:
    my $fasta_reader = new Fasta_reader($fusion_contigs_fasta);
    my %contigs_fa = $fasta_reader->retrieve_all_seqs_hash();

    ## GTF
    my %contig_to_gtf;
    {
        open(my $fh, $fusion_contigs_gtf) or die "Error, cannot open file: $fusion_contigs_gtf";
        while(<$fh>) {
            my $line = $_;
            unless ($line =~ /\w/) { next; }
            if ($line =~ /^\#/) { next; }
            my @x = split("\t", $line);
            my $contig = $x[0];
            print "gtf contig: $contig\n";
            $contig_to_gtf{$contig} .= $line;
        }
        close $fh;
    }
    


    print STDERR "-prepping individual fusion input files.\n";
    
    my %fusion_to_files;
    
    {

        ## write outputs.
        foreach my $fusion (keys(%fusion_name_to_rows)) {
            
            if (%FUSIONS_RESTRICT && ! exists $FUSIONS_RESTRICT{$fusion} ) {
                next;
            }
            
            my $fusion_output_dir = "$output_directory/$fusion";
            if (! -d $fusion_output_dir) {
                mkdir $fusion_output_dir or die "Error, cannot mkdir $fusion_output_dir";
            }
            $fusion_to_files{$fusion}->{output_dir} = $fusion_output_dir;
            
            print STDERR "-processing fusion: $fusion\n";

            eval {

                my $fusion_predictions_file = "$fusion_output_dir/fusions.tsv";
                {
                    open(my $ofh, ">$fusion_predictions_file") or die "Error, cannot write to file: $fusion_predictions_file";
                    my $delim_writer = new DelimParser::Writer($ofh, "\t", \@column_headers);
                    my @fusion_rows = @{$fusion_name_to_rows{$fusion}};
                    foreach my $fusion_row (@fusion_rows) {
                        $delim_writer->write_row($fusion_row);
                    }
                    close $ofh;
                    
                    $fusion_to_files{$fusion}->{preds_file} = $fusion_predictions_file;
                }
                
                my $fusion_genome_fasta_file = "$fusion_output_dir/fusions.fa";
                {
                    open(my $ofh, ">$fusion_genome_fasta_file") or die "Error, cannot write to $fusion_genome_fasta_file";
                    my $sequence = $contigs_fa{$fusion} or die "Error, no contig for fusion: $fusion";
                    print $ofh ">$fusion\n$sequence\n";
                    close $ofh;
                    
                    $fusion_to_files{$fusion}->{fasta_file} = $fusion_genome_fasta_file;
                    
                }
                
                my $fusion_gtf_file = "$fusion_output_dir/fusions.gtf";
                {
                    open(my $ofh, ">$fusion_gtf_file") or die "Error, cannot write to file: $fusion_gtf_file";
                    my $fusion_gtf = $contig_to_gtf{$fusion} or die "Error, no gtf info for fusion: $fusion ";;
                    print $ofh $fusion_gtf;
                    close $ofh;
                    
                    $fusion_to_files{$fusion}->{gtf_file} = $fusion_gtf_file;
                    
                }
            };
            
            if ($@) {
                print STDERR " - error writing files for fusion: $fusion_output_dir/ \n";
            }
        }
        
    }
    
    
    my $cmds_file = "$output_directory/fusion_analysis_cmds.txt";
    open(my $ofh, ">$cmds_file") or die "Error, cannot write to $cmds_file";
    
    # generate microhomology plots.
    foreach my $fusion (keys %fusion_to_files) {
        
        my $fusion_info_href = $fusion_to_files{$fusion};
        
        my $output_dir = $fusion_info_href->{output_dir};
        my $preds_file = $fusion_info_href->{preds_file};
        my $fasta_file = $fusion_info_href->{fasta_file};
        my $gtf_file = $fusion_info_href->{gtf_file};
        
        my $microhomology_file = "$output_dir/microH.dat";
        my $cmd = "($UTILDIR/misc/find_microhomologies_by_kmer_matches.pl --fasta $fasta_file --gtf $gtf_file > $microhomology_file)"
            . " && "
            . "($UTILDIR/misc/RT_artifact_inspector.Rscript --fusion_preds_tsv $preds_file --microhomologies_tsv $microhomology_file --plots_dir $output_dir)";

        print $ofh "$cmd\n";
    }
    close $ofh;
    

    ## execute them:
    if ($CPU > 0) {
        my $cmd = "ParaFly -c $cmds_file -CPU $CPU -max_retry 1 -vv";
        &process_cmd($cmd);
    }
    

    exit(0);
}


        


