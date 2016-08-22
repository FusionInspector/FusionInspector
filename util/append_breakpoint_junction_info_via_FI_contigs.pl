#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;

use lib "$FindBin::Bin/../PerlLib";
use SeqUtil;
use Fasta_reader;
use DelimParser;
use Nuc_translator;

my $usage = "usage: $0 FI_fusion_prediction_summary.dat FI.contigs.fa\n\n";

my $fusion_dat_file = $ARGV[0] or die $usage;
my $fi_contigs_fa_file = $ARGV[1] or die $usage;

my $ANCHOR_SEQ_LENGTH = 15;

main : {

    
    open (my $fh, $fusion_dat_file) or die $!;
    my $tab_reader = new DelimParser::Reader($fh, "\t");

    my @column_headers = $tab_reader->get_column_headers();

    push (@column_headers, "LeftBreakDinuc", "LeftBreakEntropy", "RightBreakDinuc", "RightBreakEntropy");

    my $tab_writer = new DelimParser::Writer(*STDOUT, "\t", \@column_headers);


    my $fasta_reader = new Fasta_reader($fi_contigs_fa_file);
    
    my %seqs = $fasta_reader->retrieve_all_seqs_hash();


    while (my $row = $tab_reader->get_row()) {
    

        my $fusion_name = $row->{'#FusionName'};
        
        my $fusion_contig_seq = $seqs{$fusion_name} or die "Error, cannot find fusion contig genomic sequence for [$fusion_name] ";
        
        my $break_left = $row->{LeftLocalBreakpoint};
        my $break_right = $row->{RightLocalBreakpoint};
        
        my ($left_splice, $left_entropy) = &examine_breakpoint_seq($break_left, $fusion_contig_seq, 'left');

        my ($right_splice, $right_entropy) = &examine_breakpoint_seq($break_right, $fusion_contig_seq, 'right');
        
        $row->{LeftBreakDinuc} = $left_splice;
        $row->{LeftBreakEntropy} = sprintf("%.4f", $left_entropy);

        $row->{RightBreakDinuc} = $right_splice;
        $row->{RightBreakEntropy} = sprintf("%.4f", $right_entropy);
        
        $tab_writer->write_row($row);
    }


    exit(0);
    

}

####
sub examine_breakpoint_seq {
    my ($coord, $fusion_contig_seq, $side) = @_;

    unless ($side eq 'left' || $side eq 'right') { 
        die "Error, cannot parse side ($side) as 'left|right' ";
    }
        
    my ($dinuc, $anchor_seq, $subseq);

    if ($side eq 'left') {
        my $coord_start = $coord - $ANCHOR_SEQ_LENGTH + 1;
        
        $subseq = substr($fusion_contig_seq, $coord_start -1, $ANCHOR_SEQ_LENGTH + 2);
        $dinuc = substr($subseq, -2);
        $anchor_seq = substr($subseq, 0, $ANCHOR_SEQ_LENGTH);
        
    }
    elsif ($side eq 'right') {
        my $coord_start = $coord - 2;
        
        $subseq = substr($fusion_contig_seq, $coord_start -1, $ANCHOR_SEQ_LENGTH + 2);
        $dinuc = substr($subseq, 0, 2);
        $anchor_seq = substr($subseq, 2);
    }
    
    #print STDERR "+:$side:$subseq|$dinuc|$anchor_seq\n";
    
        
    my $anchor_seq_entropy = &SeqUtil::compute_entropy($anchor_seq);
        

    return($dinuc, $anchor_seq_entropy);
}
