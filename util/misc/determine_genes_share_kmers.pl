#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/../../PerlLib";
use Fasta_reader;

my $usage = "usage: $0 fusion_preds.txt ref_trans.fa [KMER_LENGTH=25]\n\n";

my $fusion_preds = $ARGV[0] or die $usage;
my $trans_fa = $ARGV[1] or die $usage;

my $KMER_LENGTH = $ARGV[2] || 25;

main: {

    my $fasta_reader = new Fasta_reader($trans_fa);
    

    my %gene_id_to_seqs;
    while (my $seq_obj = $fasta_reader->next()) {
        my $header = $seq_obj->get_header();
        my ($trans_id, $gene_id, $gene_sym) = split(/\s+/, $header);
        my $sequence = $seq_obj->get_sequence();

        push (@{$gene_id_to_seqs{$gene_sym}}, $sequence);
    }
    
    
    
    open(my $fh, $fusion_preds) or die "Error, cannot open file $fusion_preds";
    while (<$fh>) {
        chomp;
        if (/^\#/) { next; }
        my $fusion_name = $_;
        my ($geneA, $geneB) = split(/--/);
        
       
        my %geneA_kmers = &get_gene_kmers($geneA, \%gene_id_to_seqs);
        my %geneB_kmers = &get_gene_kmers($geneB, \%gene_id_to_seqs);
        
                        
        my $num_kmers = &count_kmers_in_common(\%geneA_kmers, \%geneB_kmers);
                    
        print "$fusion_name\t$num_kmers\n";
    }
    
    exit(0);
}

####
sub count_kmers_in_common {
    my ($kmers_i, $kmers_j) = @_;
        
    my $count = 0;
    
    foreach my $kmer (keys %$kmers_i) {
        if (exists $kmers_j->{$kmer}) {
            $count++;
        }
    }

    return($count);
}



####
sub get_gene_kmers {
    my ($gene_id, $gene_id_to_seqs_href) = @_;

    my $seqs_aref = $gene_id_to_seqs_href->{$gene_id};
    unless ($seqs_aref) {
        die "Error, no seqs stored for gene symbol: [$gene_id] ";
    }

    my @seqs = @$seqs_aref;
    

    my %kmers;
    foreach my $seq (@seqs) {
        &seq_to_kmers(\%kmers, $seq);
    }

    return(%kmers);
}


####
sub seq_to_kmers {
    my ($kmers_href, $seq) = @_;
    
    $seq = uc $seq;

    for (my $i = 0; $i <= length($seq) - $KMER_LENGTH; $i++) {
        my $kmer = substr($seq, $i, $KMER_LENGTH);
        
        $kmers_href->{$kmer} = 1;
    }
    

    return;
}
