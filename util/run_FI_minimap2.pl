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
#  --out_prefix <string>       output prefix (default: minimap2)
#  --out_dir <string>          output directory (default: current working directory)
#  --minimap2_path <string>    full path to the minimap2 program to use.
#  --prep_reference_only       build the genome index and then stop.
#  --capture_genome_alignments reports alignments to the reference genome in addition to the fusion contigs. (for debugging)
#  --max_mate_dist <int>       maximum distance between mates (and individual introns) allowed (default: $max_mate_dist)
#  --minimap2_xtra_params <string>   extra parameters to pass on to the minimap2 aligner. Be sure to embed parameters in quotes.
#
#################################################################################################


__EOUSAGE__

    ;


my ($genome, $reads);

my $CPU = 2;

my $help_flag;

my $out_prefix = "minimap2";
my $gtf_file;
my $out_dir;
my $ADV = 0;

my $minimap2_path = "minimap2";
my $patch;
my $prep_reference_only = 0;
my $capture_genome_alignments_flag = 0;
my $samples_file;
my $minimap2_xtra_params = "";

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
             'minimap2_path=s' => \$minimap2_path,
             'prep_reference_only' => \$prep_reference_only,
             'capture_genome_alignments' => \$capture_genome_alignments_flag,

             'max_mate_dist=i' => \$max_mate_dist,

             'minimap2_xtra_params=s' => \$minimap2_xtra_params,

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


my $minimap2_prog = `which $minimap2_path`;
chomp $minimap2_prog;
unless ($minimap2_prog =~ /\w/) {
    die "Error, cannot locate minimap2 program. Be sure it's in your PATH setting.  ";
}

if ($samples_file && $reads) {
    die "Error, must specify --reads or --samples_file, not both";
}


main: {

    ## ensure all full paths
    $genome = &Pipeliner::ensure_full_path($genome);
    $gtf_file = &Pipeliner::ensure_full_path($gtf_file) if $gtf_file;

    my @read_group_samples = ();

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

        open(my $fh, "$samples_file") or die "Error, cannot open file: $samples_file";
        while (<$fh>) {
            chomp;
            my ($sample_name, $left_fq, $right_fq) = split(/\t/);
            push (@read_group_samples, $sample_name);
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
    }

    if ($out_dir) {
        unless (-d $out_dir) {
            mkdir $out_dir or die "Error, cannot mkdir $out_dir";
        }
        chdir $out_dir or die "Error, cannot cd to $out_dir";
    }


    my $combined_genome = $genome;
    my $combined_gtf = $gtf_file;

    # Handle genome patching: concatenate reference genome + fusion contigs
    if ($patch) {
        $combined_genome = "genome_w_fusion_contigs.fa";
        my $cmd = "cat $genome $patch > $combined_genome";
        &process_cmd($cmd);

        # Also combine GTF annotations if provided
        if ($gtf_file) {
            # Find patch GTF (should be in same dir as patch fasta with .gtf extension)
            my $patch_gtf = $patch;
            $patch_gtf =~ s/\.fa(sta)?$/.gtf/;

            if (-e $patch_gtf) {
                $combined_gtf = "annots_w_fusion_contigs.gtf";
                $cmd = "cat $gtf_file $patch_gtf > $combined_gtf";
                &process_cmd($cmd);
            }
        }
    }


    # Build minimap2 index
    my $mm2_index = "$combined_genome.mmi";

    unless (-e "$mm2_index.build.ok") {

        my $cmd = "$minimap2_prog -d $mm2_index $combined_genome";
        &process_cmd($cmd);

        &process_cmd("touch $mm2_index.build.ok");

        if ($prep_reference_only) {
            print STDERR "done building genome index.  stopping here.\n";
            exit(0);
        }
    }

    # Convert GTF to junction BED for minimap2
    my $splice_bed = "";
    if ($combined_gtf) {
        $splice_bed = "$combined_gtf.splice.bed";

        unless (-e "$splice_bed.ok") {
            # Use paftools.js bundled in util/ alongside this script
            my $paftools = "k8 $FindBin::RealBin/paftools.js";
            my $cmd = "$paftools gff2bed $combined_gtf > $splice_bed";
            &process_cmd($cmd);
            &process_cmd("touch $splice_bed.ok");
        }
    }


    ## run minimap2

    my $pipeliner = new Pipeliner(-verbose => 2);

    # Parse read files
    my @read_files = split(/\s+/, $reads);
    my $left_fq = $read_files[0] || "";
    my $right_fq = $read_files[1] || "";

    # Build minimap2 command
    # Note: Allow secondary alignments (default) so reads can align to both genome and fusion contigs
    # The -N parameter controls max secondary alignments (default: 5)
    my $cmd = "$minimap2_prog -ax splice -t $CPU -G $max_mate_dist -uf --MD -L ";

    if ($splice_bed && -s $splice_bed) {
        $cmd .= " --junc-bed $splice_bed ";
    }

    if ($minimap2_xtra_params) {
        $cmd .= " $minimap2_xtra_params ";
    }

    $cmd .= " $mm2_index $left_fq ";

    if ($right_fq) {
        $cmd .= " $right_fq ";
    }

    # Pipe to samtools for BAM conversion and sorting
    my $output_bam = "$out_prefix.sortedByCoord.out.bam";
    my $temp_bam = "Aligned.sortedByCoord.out.bam";

    # Filter to only fusion contig alignments before sorting (unless debugging)
    # Fusion contigs contain "--" in their reference names
    # This dramatically reduces sorting time since genome >> fusion contigs
    if ($patch && !$capture_genome_alignments_flag) {
        $cmd .= " | awk 'BEGIN{OFS=\"\\t\"} /^\@/ {print; next} \$3 ~ /--/ {print}' ";
    }

    $cmd .= " | samtools sort -@ $CPU -m 4G -o $temp_bam - ";

    $pipeliner->add_commands( new Command($cmd, "minimap2_align.ok") );


    # Handle read groups if samples file was provided
    if (@read_group_samples) {
        # Need to add read groups to the BAM file
        my $rg_added_bam = "$temp_bam.rg.bam";

        # Build read group string
        my @rg_lines;
        foreach my $sample (@read_group_samples) {
            push @rg_lines, "-r '\@RG\\tID:$sample\\tSM:$sample'";
        }
        my $rg_string = join(" ", @rg_lines);

        $cmd = "samtools addreplacerg $rg_string -o $rg_added_bam $temp_bam";
        $pipeliner->add_commands( new Command($cmd, "add_read_groups.ok") );

        $cmd = "mv $rg_added_bam $temp_bam";
        $pipeliner->add_commands( new Command($cmd, "mv_rg_bam.ok") );
    }


    # Note: Fusion contig filtering now happens in the alignment pipeline (before sorting)
    # This is much more efficient than post-alignment filtering


    # Rename to final output name
    $pipeliner->add_commands( new Command("mv $temp_bam $output_bam", "$output_bam.ok") );

    # Index the BAM file
    $pipeliner->add_commands( new Command("samtools index $output_bam", "$output_bam.bai.ok") );


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



