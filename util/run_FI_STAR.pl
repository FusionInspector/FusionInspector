#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib("$FindBin::RealBin/../PerlLib");
use Pipeliner;
use File::Basename;
use Cwd;

use Carp;
use Getopt::Long qw(:config no_ignore_case bundling pass_through);


my $max_mate_dist = 1e5;


my $usage = <<__EOUSAGE__;

#################################################################################################
#
#  Required:
#  --genome <string>           target genome to align to
#
#     --reads  <string>           fastq files. If pairs, indicate both in quotes, ie. "left.fq right.fq"
#     or
#     --samples_file <string>     samples file (format:  sample(tab)/path/left.fq(tab)/path/right.fq)
#
#  --patch <string>            genomic targets to patch the genome fasta with.
#
#
#  Optional:
#  -G <string>                 GTF file for incorporating reference splice site info.
#  --CPU <int>                 number of threads (default: 2)
#  --out_prefix <string>       output prefix (default: star)
#  --out_dir <string>          output directory (default: current working directory)
#  --star_path <string>        full path to the STAR program to use.
#  --prep_reference_only       build the genome index and then stop.
#  --only_fusion_reads         restrict alignments only to the fusion-supporting reads
#  --capture_genome_alignments reports alignments to the reference genome in addition to the fusion contigs. (for debugging)
#  --chim_search               include Chimeric.junction outputs
#  --max_mate_dist <int>       maximum distance between mates (and individual introns) allowed (default: $max_mate_dist)
#  --no_splice_score_boost     do not augment alignment score for spliced alignments 
#  --STAR_xtra_params <string>   extra parameters to pass on to the STAR aligner. Be sure to embed parameters in quotes.
#
#################################################################################################


__EOUSAGE__

    ;


my ($genome, $reads);

my $CPU = 2;

my $help_flag;

my $out_prefix = "star";
my $gtf_file;
my $out_dir;
my $ADV = 0;

my $star_path = "STAR";
my $patch;
my $prep_reference_only = 0;
my $only_fusion_reads_flag = 0;
my $capture_genome_alignments_flag = 0;
my $chim_search;
my $samples_file;
my $no_splice_score_boost = 0;
my $STAR_xtra_params = "";

&GetOptions( 'h' => \$help_flag,
             'genome=s' => \$genome,
             'patch=s' => \$patch,
             
             'reads=s' => \$reads,
             'samples_file=s' => \$samples_file,
             
             'CPU=i' => \$CPU,
             'out_prefix=s' => \$out_prefix,
             'G=s' => \$gtf_file,
             'out_dir=s' => \$out_dir,
             'ADV' => \$ADV,
             'star_path=s' => \$star_path,
             'prep_reference_only' => \$prep_reference_only,
             'only_fusion_reads' => \$only_fusion_reads_flag,
             'capture_genome_alignments' => \$capture_genome_alignments_flag,
             'chim_search' => \$chim_search,

             'max_mate_dist=i' => \$max_mate_dist,

             'no_splice_score_boost' => \$no_splice_score_boost,

             'STAR_xtra_params=s' => \$STAR_xtra_params,

    );


unless ($genome && ($reads || $samples_file) ) {
    die $usage;
}

if ($help_flag) {
    die $usage;
}

if (@ARGV) {
    die "Error, cannot recognize opts: @ARGV";
}


my $star_prog = `which $star_path`;
chomp $star_prog;
unless ($star_prog =~ /\w/) {
    die "Error, cannot locate STAR program. Be sure it's in your PATH setting.  ";
}

if ($samples_file && $reads) {
    die "Error, must specify --reads or --samples_file, not both";
}


main: {
	
    ## ensure all full paths
    $genome = &Pipeliner::ensure_full_path($genome);
    $gtf_file = &Pipeliner::ensure_full_path($gtf_file) if $gtf_file;

    my $read_group_ids = "";
    
    if ($reads) {
        my @read_files = split(/\s+/, $reads);
        foreach my $read_file (@read_files) {
            if ($read_file) {
                $read_file = &Pipeliner::ensure_full_path($read_file);
            }
        }
        $reads = join(" ", @read_files);
    }
    else {

        my @left_fqs;
        my @right_fqs;
        my @read_groups;
        
        open(my $fh, "$samples_file") or die "Error, cannot open file: $samples_file";
        while (<$fh>) {
            chomp;
            my ($sample_name, $left_fq, $right_fq) = split(/\t/);
            push (@read_groups, "ID:$sample_name");
            unless (-s $left_fq) {
                confess "Error, cannot locate file: $left_fq";
            }
            push (@left_fqs, &Pipeliner::ensure_full_path($left_fq));

            if ($right_fq) {
                unless (-s $right_fq) {
                    confess "Error, cannot locate file: $right_fq";
                }
                push (@right_fqs, &Pipeliner::ensure_full_path($right_fq));
            }
        }
        close $fh;

        $reads = join(",", @left_fqs);
        if (@right_fqs) {
            $reads .= " " . join(",", @right_fqs);
        }
        $read_group_ids = join(" , ", @read_groups);
    }
    
    if ($out_dir) {
        unless (-d $out_dir) {
            mkdir $out_dir or die "Error, cannot mkdir $out_dir";
        }
        chdir $out_dir or die "Error, cannot cd to $out_dir";
    }

    my $MIN_RAM = 1024**3; # 1G
    my $genome_size = -s $genome;
    my $estimated_ram = $genome_size * 15;
    
    if ($estimated_ram < $MIN_RAM) {
        $estimated_ram = $MIN_RAM;
    }
    
    my $star_index = "$genome.star.idx";
    if (! -d $star_index) {
        
        ## build star index
        unless (-d $star_index) {
            mkdir($star_index) or die "Error, cannot mkdir $star_index";
        }

        # from Alex D.:
        # scale down the --genomeSAindexNbases parameter as log2(GenomeLength)/2 - 1
        
    
        my $genomeSAindexNbases = int(log($genome_size) / log(2) / 2); # close enough.
    

    
        my $cmd = "$star_prog --runThreadN $CPU --runMode genomeGenerate --genomeDir $star_index "
            . " --genomeFastaFiles $genome "
            . " --genomeSAindexNbases $genomeSAindexNbases"
            . " --limitGenomeGenerateRAM $estimated_ram "; #40419136213 ";
        if ($gtf_file) {
            $cmd .= " --sjdbGTFfile $gtf_file "
                . " --sjdbOverhang 150 ";
            
        }
    
        &process_cmd($cmd);
        
        &process_cmd("touch $star_index/build.ok");
        
        if ($prep_reference_only) {
            print STDERR "done building genome index.  stopping here.\n";
            exit(0);
        }
    }

        
    ## run STAR
        my @tmpfiles;
    
    my $pipeliner = new Pipeliner(-verbose => 2);
    
    my $cmd = "$star_prog "
        . " --runThreadN $CPU "
        . " --genomeDir $star_index "
        . " --outSAMtype BAM SortedByCoordinate "
        . " --twopassMode None "
        . " --alignSJDBoverhangMin 10 "
        . " --genomeSuffixLengthMax 10000"
        . " --limitBAMsortRAM $estimated_ram "  #20000000000";
        . " --alignInsertionFlush Right  "
        . " --alignMatesGapMax $max_mate_dist "
        . " --alignIntronMax $max_mate_dist "
        . " --outSAMattributes All ";
    
    if (length($reads) > 1000) {
        my $star_params_file = "star_params.$$.txt";
        open(my $ofh, ">$star_params_file") or die "Error, cannot write to file: $star_params_file";
        print $ofh "readFilesIn $reads\n";
        if ($read_group_ids) {
            print $ofh "outSAMattrRGline $read_group_ids\n";
        }
        close $ofh;

        $cmd .= " --parametersFiles $star_params_file ";
    }
    else {
        $cmd .= " --readFilesIn $reads ";
        if ($read_group_ids) {
            $cmd .= " --outSAMattrRGline $read_group_ids ";
        }
    }
    

    if ($chim_search) {
        $cmd .= " --chimJunctionOverhangMin 12 "
            .  " --chimSegmentMin 12 "
            .  " --chimSegmentReadGapMax 3 "
    }
    
    if ($patch) {
        
        $cmd .= " --genomeFastaFiles $patch ";
        
        if ($only_fusion_reads_flag ) {
            $cmd .= " --outSAMfilter KeepOnlyAddedReferences ";
        }
        elsif ($capture_genome_alignments_flag) {
            # no op
            # by default, reports all alignments
        }
        else {
            ## ** the FusionInspector default setting ** 
            # retains all alignments to the contigs, not just the fusion evidence.
            $cmd .= " --outSAMfilter KeepAllAddedReferences ";
        }
    }
    
    if ($gtf_file) {
        $cmd .= " --sjdbGTFfile $gtf_file ";
    }

    if ($no_splice_score_boost) {
        $cmd .= " --alignSJstitchMismatchNmax 5 5 5 5 "; # same for all breakpoints
        $cmd .= " --scoreGapGCAG -8 "; # match the other dinuc pair penalties.
        
    } else {
        $cmd .= " --alignSJstitchMismatchNmax 5 -1 5 5 ";  #which allows for up to 5 mismatches for non-canonical GC/AG, and AT/AC junctions, and any number of mismatches for canonical junctions (the default values 0 -1 0 0 replicate the old behavior (from AlexD)
        
        $cmd .= " --scoreGapNoncan -6 "; # default is -8 ... want to penalize slightly less, othwerise can get mismatches near false splice junctions in otherwise optimal alignments.
    }
    
    if ($reads =~ /\.gz$/) {
        $cmd .= " --readFilesCommand 'gunzip -c' ";
    }
    
    if ($STAR_xtra_params) {
        $cmd .= " $STAR_xtra_params ";
    }
    
    $pipeliner->add_commands( new Command($cmd, "star_align.ok") );
    
    
    my $bam_outfile = "Aligned.sortedByCoord.out.bam";
    my $renamed_bam_outfile = "$out_prefix.sortedByCoord.out.bam";
    $pipeliner->add_commands( new Command("mv $bam_outfile $renamed_bam_outfile", "$renamed_bam_outfile.ok") );
    
    
    $pipeliner->add_commands( new Command("samtools index $renamed_bam_outfile", "$renamed_bam_outfile.bai.ok") );
    
    
    $pipeliner->run();

    
	exit(0);
}



####
sub process_cmd {
	my ($cmd) = @_;
	
	print STDERR "CMD: $cmd\n";
	#return;

	my $ret = system($cmd);
	if ($ret) {
		die "Error, cmd: $cmd died with ret ($ret)";
	}

	return;
}



