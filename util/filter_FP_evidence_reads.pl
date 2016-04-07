#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Long qw(:config no_ignore_case bundling pass_through);

use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Data::Dumper;
use Pipeliner;
use Fastq_reader;
use SAM_reader;
use SAM_entry;
use Cwd;


unless ($ENV{TRINITY_HOME}) {
    die "Error, must have Trinity installed and env var TRINITY_HOME set to the Trinity installation directory";
}

my $CPU = 1;

my $usage = <<__EOUSAGE__;

###############################################################################################
#
#  --left_fq <string>             left fastq filename
#  
#  --right_fq <string>            right fastq filename
#
#  --cdna_fa <string>             reference cdna fasta filename
#
#  --fusion_summary <string>      fusion summary file (pre-filtering)
#
#  --genome_lib_dir <string>      genome lib dir
#
#  *Optional
#
#  --tmpdir <string>              dir for tmpfiles (default is curr dir)
#
#  --CPU <int>                    multithreading for trimmomatic (default: $CPU)
#
#  --quality_trim                 quality trim the evidence reads
#
###############################################################################################


__EOUSAGE__

    ;

my $left_fq;
my $right_fq;
my $cdna_fa;
my $fusion_summary;
my $help_flag;
my $tmpdir = cwd();
my $DO_QUALITY_TRIMMING;

my $MAX_KMER_ABUNDANCE = 1;

my $genome_lib_dir;

&GetOptions ( 'h' => \$help_flag,
              
              'left_fq=s' => \$left_fq,
              'right_fq=s' => \$right_fq,
              'cdna_fa=s' => \$cdna_fa,
              'fusion_summary=s' => \$fusion_summary,
    
              'tmpdir=s' => \$tmpdir,
              'CPU=i' => \$CPU,
       
              'quality_trim' => \$DO_QUALITY_TRIMMING,
              
              'genome_lib_dir=s' => \$genome_lib_dir,

    );


if ($help_flag) {
    die $usage;
}

unless ($left_fq && $right_fq && $cdna_fa && $fusion_summary && $genome_lib_dir) {
    die $usage;
}

my $TRINITY_HOME = $ENV{TRINITY_HOME} or die "Error, need env var TRINITY_HOME set to Trinity installation directory";
my $TRIMMOMATIC_JAR = "$TRINITY_HOME/trinity-plugins/Trimmomatic/trimmomatic.jar";
unless (-s $TRIMMOMATIC_JAR) {
    die "Error, cannot locate trimmomatic jar file at $TRIMMOMATIC_JAR";
}
my $TRIMMOMATIC_SETTINGS = "SLIDINGWINDOW:4:5 LEADING:5 TRAILING:5 MINLEN:25";

main: {


    print STDERR "-parsing junction and spanning reads info.\n";
    my %junction_reads;
    my %spanning_frags;
    &parse_junction_and_spanning_reads($fusion_summary, \%junction_reads, \%spanning_frags);

    print STDERR "-identified " . scalar(keys %junction_reads) . " junction reads.\n";
    print STDERR "-identified " . scalar(keys %spanning_frags) . " spanning frags.\n";
    

    
    ##
    ## Examine/filter the junction reads
    ##
    
    print STDERR "-extracting junction reads from fq file\n";
    
    my $junc_reads_fq = "$tmpdir/tmp.junc_reads.fq";
    ## examine junction reads
    
    ## get the reads
    
    if ($DO_QUALITY_TRIMMING) {
        
        my $pre_Qtrimmed_reads = "$junc_reads_fq.pre_Qtrimmed";
        
        &extract_junction_reads(\%junction_reads, $left_fq, $right_fq, $pre_Qtrimmed_reads);
        
        &quality_trim_junc_reads(\%junction_reads, $pre_Qtrimmed_reads, $junc_reads_fq);
    }
    else {
        &extract_junction_reads(\%junction_reads, $left_fq, $right_fq, $junc_reads_fq);
    }
    
    
    ## align the reads to the cdna fasta file
    my $junc_reads_sam = "$junc_reads_fq.sam";
    my $sam_chkpt = "$junc_reads_sam.ok";
    my $bowtie_cmd = "bowtie -q -S --sam-nohead --best --chunkmbs 512 $cdna_fa  $junc_reads_fq > $junc_reads_sam";
    &process_cmd($bowtie_cmd) unless (-e $sam_chkpt);
    
    &process_cmd("touch $sam_chkpt");
    
    ## eliminate those junction reads that do actually align ok to reference transcripts
    my $number_reads_mapped = 0;
    my $sam_reader = new SAM_reader($junc_reads_sam);
    while (my $sam_entry = $sam_reader->get_next()) {
        if (! $sam_entry->is_query_unmapped()) {
            my $full_read_name = $sam_entry->get_read_name();
            $number_reads_mapped++;
            
            if (delete $junction_reads{$full_read_name}) {
                print STDERR "-deleted junction read: $full_read_name\n";
            }
            else {
                print STDERR Dumper(\%junction_reads);
                die "Error, cannot delete junction read: $full_read_name\n";
            }
            
            
        }
    }
    print STDERR "- $number_reads_mapped junction reads excluded due to mapping ok to reference transcripts\n";
    
    
    
    ##
    ## Examine/filter the spanning frags
    ##
    
    print STDERR "-extracting spanning frags from fq files\n";
    
    my $span_reads_left_fq = "$tmpdir/tmp.span_reads.left.fq";
    my $span_reads_right_fq = "$tmpdir/tmp.span_reads.right.fq";
    
    
    if ($DO_QUALITY_TRIMMING) {
        
        my $pre_Qtrimmed_left_fq = "$span_reads_left_fq.preQtrimmed.fq";
        my $pre_Qtrimmed_right_fq = "$span_reads_right_fq.preQtrimmed.fq";
        
        &extract_spanning_reads(\%spanning_frags, $left_fq, $pre_Qtrimmed_left_fq);
        &extract_spanning_reads(\%spanning_frags, $right_fq, $pre_Qtrimmed_right_fq);
        
        &quality_trim_spanning_reads(\%spanning_frags, 
                                     $pre_Qtrimmed_left_fq, $pre_Qtrimmed_right_fq, 
                                     $span_reads_left_fq, $span_reads_right_fq);
    }
    else {
        
        &extract_spanning_reads(\%spanning_frags, $left_fq, $span_reads_left_fq);
        &extract_spanning_reads(\%spanning_frags, $right_fq, $span_reads_right_fq);
        
    }
    
    ## align reads, identify those that actually align as proper pairs.
    my $span_frags_sam = "$tmpdir/tmp.span_reads.sam";
    my $sam_frags_chkpt = "$span_frags_sam.ok";
    
    my $bowtie_span_cmd = "bowtie -q -S --sam-nohead --best --chunkmbs 512 $cdna_fa -1 $span_reads_left_fq -2 $span_reads_right_fq > $span_frags_sam";
    &process_cmd($bowtie_span_cmd) unless (-e $sam_frags_chkpt);
    
    ## eliminate those spanning frags that do actually align ok to reference transcripts
    my $number_frags_mapped = 0;
    my $sam_frag_reader = new SAM_reader($span_frags_sam);
    while (my $sam_entry = $sam_frag_reader->get_next()) {
        if (! $sam_entry->is_query_unmapped()) {
            my $core_read_name = $sam_entry->get_core_read_name();
            if ( (exists $spanning_frags{$core_read_name}) && $sam_entry->is_proper_pair() ) {
                
                $number_frags_mapped++;
                delete $spanning_frags{$core_read_name} or die "Error, could not delete $core_read_name from spanning_frags";
            }
        }
    }
    &process_cmd("touch $sam_frags_chkpt");
    
    print STDERR "- $number_frags_mapped spanning frags excluded due to mapping ok to reference transcripts\n";
    
    


    #########################################
    ## Filter any reads that are repetitive
    #########################################
    
    ## TODO:  this needs to be parameterized, hardcoding now for testing purposes.
    
    &filter_repetitive_reads(\%junction_reads, \%spanning_frags, $junc_reads_fq, $span_reads_left_fq, $span_reads_right_fq);
    

    
    ##
    ## report the adjusted fusion summary, removing the false evidence:
    ##
    
    print STDERR "-post-filtering, have " . scalar(keys %junction_reads) . " junction reads.\n";
    print STDERR "-post-filtering, have " . scalar(keys %spanning_frags) . " spanning frags.\n";
    
    &exclude_FP_junction_and_spanning_reads($fusion_summary, \%junction_reads, \%spanning_frags);
    

    exit(0);

}

####
sub process_cmd {
    my ($cmd) = @_;

    print STDERR "CMD: $cmd\n";
    my $ret = system($cmd);
    if ($ret) {
        confess "Error, cmd: $cmd died with ret $ret";
    }
    return;
}



####
sub filter_repetitive_reads {
    my ($junction_reads_href, $spanning_frags_href, 
        @fq_files) = @_;

    my $kmer_repeat_info = "$tmpdir/repeat_kmer_info.dat";

    my $filter_repetitive_checkpoint = "$tmpdir/filter_repetitive.ok";
    unless (-e $filter_repetitive_checkpoint) {

        
        my $tmp_fa = "$tmpdir/tmp.reads.fa";
        if (-e $tmp_fa) {
            unlink($tmp_fa) or die "error, cannot remove file $tmp_fa";
        }
        

        foreach my $fq_file (@fq_files) {
            my $cmd = "$ENV{TRINITY_HOME}/util/support_scripts/fastQ_to_fastA.pl -I $fq_file >> $tmp_fa";
            &process_cmd($cmd);
        }
        
        my $kmers_file = "$genome_lib_dir/ref_genome.jf.min2.kmers";
        my $cmd = "$ENV{TRINITY_HOME}/Inchworm/bin/fastaToKmerCoverageStats --kmers $kmers_file --reads $tmp_fa --DS --capture_coverage_info > $kmer_repeat_info";
        &process_cmd($cmd);
    }
    
    my $count_removed = 0;
    
    open (my $fh, $kmer_repeat_info) or die "Error, cannot open file $kmer_repeat_info";
    while (<$fh>) {
        chomp;
        my @x = split(/\t/);
        my $median_kmer_abundance = $x[0];
        my $read_name = $x[4];
        
        if ($median_kmer_abundance > $MAX_KMER_ABUNDANCE) {
            
            $count_removed++;

            ## exclude
            if (exists $junction_reads_href->{$read_name}) {
                delete $junction_reads_href->{$read_name};
            }
            
            $read_name =~ s|/[12]$||;
            if (exists $spanning_frags_href->{$read_name}) {
                delete $spanning_frags_href->{$read_name};
            }
        }
    }

    print STDERR "$count_removed reads removed as repetitive.\n";
    
    return;
}




####
sub extract_junction_reads {
    my ($junction_reads_href, $left_fq, $right_fq, $junc_reads_fq) = @_;

    my %junc_reads = %$junction_reads_href;

    my $ofh;
    my @fq_files;

    my $checkpoint = "$junc_reads_fq.ok";

    ## re-use earlier generated output file if available.
    if (-s $junc_reads_fq && -e $checkpoint) {
        @fq_files = ($junc_reads_fq);
        print STDERR "-resuming using $junc_reads_fq file\n";
    }
    else {
        @fq_files = ($left_fq, $right_fq);
        open ($ofh, ">$junc_reads_fq") or confess "Error, cannot write to file $junc_reads_fq";
    }
    
    foreach my $fq_file (@fq_files) {
        
        print STDERR "\n\n-parsing $fq_file\n";
        my $counter = 0;
        
        my $fastq_reader = new Fastq_reader($fq_file);
        while (my $fq_obj = $fastq_reader->next()) {
            my $full_read_name = $fq_obj->get_full_read_name();
            if (exists $junc_reads{$full_read_name}) {
                my $fq_record = $fq_obj->get_fastq_record();
                
                print $ofh $fq_record if $ofh;
                delete $junc_reads{$full_read_name};
                
            }
            $counter++;
            if ($counter % 10000 == 0) {
                print STDERR "\r[$counter]  ";
            }
        }
    }
    print STDERR "\n\n";
    
    if ($ofh) {
        close $ofh;
        &process_cmd("touch $checkpoint");
    }
    
    if (%junc_reads) {
        confess "Error, didn't retrieve all junction reads.  Missing: " . Dumper(\%junc_reads);
    }

    return;
}

####
sub extract_spanning_reads {
    my ($spanning_frags_href, $fastq_input, $fastq_output) = @_;

    my %frags = %$spanning_frags_href;
    
    my $ofh;

    # reuse existing output file if available.
    my $fastq_file_to_read;
    
    my $checkpoint = "$fastq_output.ok";
    if (-s $fastq_output && -e $checkpoint) {
        $fastq_file_to_read = $fastq_output;
        print STDERR "-resuming using $fastq_output\n";
    }
    else {
        $fastq_file_to_read = $fastq_input;
        open ($ofh, ">$fastq_output") or die "error, cannot write to file $fastq_output";        
    }

    print STDERR "\n-reading $fastq_file_to_read\n\n";
    my $fastq_reader = new Fastq_reader($fastq_file_to_read);
    
    while (my $fq_obj = $fastq_reader->next()) {
        my $counter = 0;
        my $core_read_name = $fq_obj->get_core_read_name();
        if (exists $frags{$core_read_name}) {
            if ($ofh) {
                my $fq_record = $fq_obj->get_fastq_record();
                print $ofh $fq_record;
            }
            delete($frags{$core_read_name});
        }
        $counter++;
        if ($counter % 10000 == 0) {
            print "\r[$counter]  ";
        }
    }
    print STDERR "\n\n";
    

    if ($ofh) {
        close $ofh;
        &process_cmd("touch $checkpoint");
    }
    
    if (%frags) {
        confess "Error, missing reads from file $fastq_input: " . Dumper(\%frags);
    }
    
    return;
}


####
sub parse_junction_and_spanning_reads {
    my ($fusion_summary, $junction_reads_href, $spanning_frags_href) = @_;

    open (my $fh, $fusion_summary) or die "Error, cannot open file $fusion_summary";
    while (<$fh>) {
        chomp;
        if (/^\#/) { next; }
        my @x = split(/\t/);
        my $junc_reads_list = $x[10];
        my $spanning_frag_list = $x[11];
        
        foreach my $junc_read (split(/,/, $junc_reads_list)) {
            $junction_reads_href->{$junc_read}++;
        }
        foreach my $spanning_frag (split(/,/, $spanning_frag_list)) {
            $spanning_frags_href->{$spanning_frag}++;
        }
    }
    close $fh;

    return;
}


####
sub exclude_FP_junction_and_spanning_reads {
    my ($fusion_summary, $junction_reads_href, $spanning_frags_href) = @_;
    
    open (my $ofh, ">$fusion_summary.filt") or die $!;
    
    open (my $fh, $fusion_summary) or die "Error, cannot open file $fusion_summary";
    while (<$fh>) {
        if (/^\#/) { print; next; }
        chomp;
        my @x = split(/\t/);
        my $junc_reads_list = $x[10];
        my $spanning_frag_list = $x[11];

        my @adj_junc_reads;
        foreach my $junc_read (split(/,/, $junc_reads_list)) {
            if (exists $junction_reads_href->{$junc_read}) {
                push (@adj_junc_reads, $junc_read);
            }
        }

        my @adj_spanning_frags;
        foreach my $spanning_frag (split(/,/, $spanning_frag_list)) {
            if (exists $spanning_frags_href->{$spanning_frag}) {
                push (@adj_spanning_frags, $spanning_frag);
            }
        }

        ## adjust the stats based on the adjusted evidence counts.
        
        $x[10] = join(",", @adj_junc_reads);
        $x[11] = join(",", @adj_spanning_frags);
        
        my $orig_junc_read_count = $x[7];
        my $orig_span_frag_count = $x[8];
        
        my $num_junction_reads = $x[7] = scalar(@adj_junc_reads);
        my $num_spanning_reads = $x[8] = scalar(@adj_spanning_frags);
        
        my $num_left_contrary_reads = $x[12];
        my $num_right_contrary_reads = $x[14];

        my $TAF_left = ($num_junction_reads + $num_spanning_reads) 
            / 
            ($num_left_contrary_reads $num_junction_reads + $num_spanning_reads);
        
        $TAF_left = sprintf("%.2f", $TAF_left);
        
        my $TAF_right = ($num_junction_reads + $num_spanning_reads) 
            / 
            ($num_right_contrary_reads + $num_junction_reads + $num_spanning_reads);
        
        $TAF_right = sprintf("%.2f", $TAF_right);
        
        $x[16] = $TAF_left;
        $x[17] = $TAF_right;
        
        
        my $pct_filtered_junction = 0;
        if ($orig_junc_read_count > 0) {
            $pct_filtered_junction = sprintf("%.2f", ($orig_junc_read_count - $num_junction_reads) / $orig_junc_read_count * 100);
        }
        
        my $pct_filtered_spanning = 0;
        if ($orig_span_frag_count > 0) {
            $pct_filtered_spanning = sprintf("%.2f", ($orig_span_frag_count - $num_spanning_reads) / $orig_span_frag_count * 100);
        }
        

        if ($num_junction_reads == 0 && $x[6] ne "NO_JUNCTION_READS_IDENTIFIED") {
            # reset it
            $x[6] = "NO_JUNCTION_READS_IDENTIFIED";
        }
        
        #add % filtered annotation
        if ($pct_filtered_junction > 0 || $pct_filtered_spanning > 0) {
            # add to annotations
            if ($x[18] eq ".") {
                $x[18] = "";
            }
            else {
                $x[18] .= ",";
            }
            
            $x[18] .= "PctFiltJ[$pct_filtered_junction],PctFiltS[$pct_filtered_spanning]";
        }
        
        # report entry
        my $outline = join("\t", @x) . "\n"; 
        
        if ($num_junction_reads > 0 || $num_spanning_reads > 0) {
            print $outline;
            print $ofh $outline;
        }
        else {
            print $ofh "#" . $outline;
        }
        
    }
    close $fh;
    close $ofh;

    return;
}


####
sub quality_trim_junc_reads {
    my ($junction_reads_href, $pre_Qtrimmed_reads, $junc_reads_fq) = @_;
    
    my $cmd = "java -Xmx1G -jar $TRIMMOMATIC_JAR SE -threads $CPU -phred33 $pre_Qtrimmed_reads $junc_reads_fq $TRIMMOMATIC_SETTINGS";
    &process_cmd($cmd);

    my %trimmed_reads;
    my $fastq_reader = new Fastq_reader($junc_reads_fq) or confess "Error, cannot open file $junc_reads_fq";
    while (my $fq_obj = $fastq_reader->next()) {

        my $read_name = $fq_obj->get_full_read_name();
        $trimmed_reads{$read_name} = 1;
    }

    my @junc_read_accs = keys %$junction_reads_href;
    foreach my $junc_read_acc (@junc_read_accs) {
        if (exists $trimmed_reads{$junc_read_acc}) {
            delete $trimmed_reads{$junc_read_acc};
            # retain as junc read
        }
        else {
            # no longer keep in junc reads
            delete $junction_reads_href->{$junc_read_acc} or confess "Error, wasnt able to delete $junc_read_acc from junction read list";
        }
    }

    if (%trimmed_reads) {
        confess "Error, some trimmed reads weren't recognized among the junction read set: " . Dumper(\%trimmed_reads);
    }

    return;
}


####
sub quality_trim_spanning_reads {
    my ($spanning_frags_href, $pre_Qtrim_left_fq, $pre_Qtrim_right_fq, $span_left_fq, $span_right_fq) = @_;

    my $cmd = "java -Xmx1G -jar $TRIMMOMATIC_JAR PE -threads $CPU -phred33 "
        . " $pre_Qtrim_left_fq $pre_Qtrim_right_fq "
        . " $span_left_fq $span_left_fq.U "
        . " $span_right_fq $span_right_fq.U "
        . " $TRIMMOMATIC_SETTINGS ";

    &process_cmd($cmd);

    my %trimmed_span_frags;
    my $fastq_reader = new Fastq_reader($span_left_fq);
    while (my $fq_obj = $fastq_reader->next() ) {

        my $read_name = $fq_obj->get_core_read_name();
        $trimmed_span_frags{$read_name} = 1;
    
    }
    
    my @span_frag_accs = keys %$spanning_frags_href;
    foreach my $span_frag_acc (@span_frag_accs) {
        if (exists $trimmed_span_frags{$span_frag_acc}) {
            delete $trimmed_span_frags{$span_frag_acc};
        }
        else {
            delete $spanning_frags_href->{$span_frag_acc} or confess "Error, wasant able to delete $span_frag_acc from spanning frag list";
        }
    }
    
    if (%trimmed_span_frags) {
        confess "Error, some trimmed spanning frags werent recognized among the spanning frag set: " . Dumper(\%trimmed_span_frags);
    }


    return;
}
