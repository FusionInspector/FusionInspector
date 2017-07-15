#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use SAM_reader;
use SAM_entry;
use Data::Dumper;
use Carp;
use DelimParser;

my $usage = "usage: $0 read_names.accs  fileA.bam,fileB.bam,...\n\n";

my $read_name_accs_list = $ARGV[0] or die $usage;
my $bam_file_listing = $ARGV[1] or die $usage;

main: {
    
    my %reads_want;
    {
        open(my $fh, $read_name_accs_list) or die "Error, cannot open file $read_name_accs_list";
        my $tab_reader = new DelimParser::Reader($fh, "\t");
        
        while (my $row = $tab_reader->get_row()) {
            my $geneA = $row->{LeftGene};
            my $geneB = $row->{RightGene};
            my $reads_list = $row->{JunctionReads};

            $geneA =~ s/\^.*$//;
            $geneB =~ s/\^.*$//;
                        
            my $fusion_contig = "$geneA--$geneB";
            foreach my $read_name (split(/,/, $reads_list)) {
                $read_name =~ s/\/[12]$//;
                $reads_want{"$fusion_contig|$read_name"} = 1;
            }
        }
                
    }

    my %seen;

    foreach my $bam_file (split(/,/, $bam_file_listing) ) {
        
        my $sam_reader = new SAM_reader($bam_file);
        while (my $sam_entry = $sam_reader->get_next()) {
            my $scaffold = $sam_entry->get_scaffold_name();
            my $read_name = $sam_entry->get_core_read_name();
            #my $read_name = $sam_entry->reconstruct_full_read_name();
            $read_name = "$scaffold|$read_name";
            if ($reads_want{$read_name}) {
                print $sam_entry->get_original_line() . "\n";
                $seen{$read_name} = 1;
            }
        }

    }
    
    foreach my $read_name (keys %seen) {
        delete $reads_want{$read_name};
    }
    
    if (%reads_want) {
        confess "Error, missing junction reads: " . Dumper(\%reads_want);
        
    }
    
    exit(0);
}


