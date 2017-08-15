#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use DelimParser;
use Data::Dumper;
use Carp;

my $usage = "\n\tusage: $0 fusion_summary_file TrinGG_fusion_gff3\n\n";

my $fusion_summary_file = $ARGV[0] or die $usage;
my $trin_fusions_gff3 = $ARGV[1] or die $usage;


main: {

    my %breakpoint_info;
    {
        open (my $fh, $trin_fusions_gff3) or die "Error, cannot open file $trin_fusions_gff3";
        while (<$fh>) {
            chomp;
            if (/^\#TrinityFusionTranscript/) {
                chomp;
                my @x = split(/\t/);
                my ($token, $trin_acc, $contig_brkpt) = @x;

                push (@{$breakpoint_info{$contig_brkpt}}, $trin_acc);
            }
        }
        close $fh;
    }
    
    open (my $fh, $fusion_summary_file) or die "Error, cannot open file $fusion_summary_file";
    my $tabreader = new DelimParser::Reader($fh, "\t");
    my @column_headers = $tabreader->get_column_headers();

    push (@column_headers, "TrinGG_Fusion");

    my $tabwriter = new DelimParser::Writer(*STDOUT, "\t", \@column_headers);

    while (my $row = $tabreader->get_row()) {
        
        my $geneA = $row->{'LeftGene'} or confess "error, no LeftGene in " . Dumper($row);
        my $local_brkpt_A = $row->{'LeftLocalBreakpoint'} or confess "error, no LeftLocalBreakpoint in " . Dumper($row);

        my $geneB = $row->{'RightGene'} or confess "error, no RightGene in " . Dumper($row);
        my $local_brkpt_B = $row->{'RightLocalBreakpoint'} or confess "error, no RightLocalBreakpoint in " . Dumper($row);
        
        $geneA =~ s/\^.*$//;
        $geneB =~ s/\^.*//;

        $row->{'TrinGG_Fusion'} = '.';
                
        my $fusion_name = "$geneA--$geneB";
        my $fusion_brkpt = "$fusion_name:$local_brkpt_A-$local_brkpt_B";
        if (my $trin_list_aref = $breakpoint_info{$fusion_brkpt}) {
            $row->{'TrinGG_Fusion'} = join(",", @$trin_list_aref);
        }
        $tabwriter->write_row($row);
    }
    close $fh;

    

    exit(0);
}


            
