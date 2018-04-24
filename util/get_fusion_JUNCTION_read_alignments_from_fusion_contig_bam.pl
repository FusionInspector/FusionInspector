#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use SAM_reader;
use SAM_entry;
use DelimParser;
use Carp;
use Data::Dumper;
use Set::IntervalTree;
use SeqUtil;
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);

my $DEBUG = 0;

my $usage = <<__EOUSAGE__;

###############################################################
#
# Required:
#
#  --junction_read_info <string>   finspector.star.cSorted.dupsMarked.bam.fusion_junction_info
#
#  --bam <string>                  read_alignments.bam
#
###############################################################


__EOUSAGE__

    ;


my $junction_read_info_file;
my $bam_file;
my $help_flag;

&GetOptions('help|h' => \$help_flag,
          
            'junction_read_info=s' => \$junction_read_info_file,
            'bam=s' => \$bam_file,
            
            'd' => \$DEBUG,
  );


if ($help_flag) {
    die $usage;
}

unless ($junction_read_info_file && $bam_file) {
    die $usage;
}


main: {

    
    my %alignments_want;
    
    {
        open(my $fh, $junction_read_info_file) or die "Error, cannot open file: $junction_read_info_file";
        my $delim_parser = new DelimParser::Reader($fh, "\t");
        while (my $row = $delim_parser->get_row()) {
            my $contig = $row->{FusionName};
            unless (defined $contig) {
                die "Error, not finding FusionName info in " . Dumper($row);
            }
            my $junction_reads = $row->{JunctionReads};
            my @reads = split(/,/, $junction_reads);
            foreach my $read (@reads) {
                $read =~ s/\/[12]$//;
                my $read_scaff_token = join("$;", $contig, $read);
                $alignments_want{$read_scaff_token} = 1;
            }
        }
        close $fh;
    }
    

    {
        # capture the alignments involving identified junction / split reads:

        my %must_find = %alignments_want;
        
        my $sam_reader = new SAM_reader($bam_file);
        my $counter = 0;
        while (my $sam_entry = $sam_reader->get_next()) {

            $counter++;

            my $contig_name = $sam_entry->get_scaffold_name();
            my $core_read_name = $sam_entry->get_core_read_name();

            my $alignment_token = join("$;", $contig_name, $core_read_name);
            
            # want both reads in the pair
            if ($alignments_want{$alignment_token}) {
                print $sam_entry->get_original_line() . "\n";
                delete($must_find{$alignment_token}) if (exists $must_find{$alignment_token});
            }
            
        }
    
    }

    exit(0);
    
}

