#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);
use lib ($ENV{EUK_MODULES});
use Overlap_piler;
use Fasta_reader;

my $KMER_SIZE = 10;

my $usage = <<__EOUSAGE__;


############################################
#
# Required:
#
#  --fasta <string>    finspector.fasta
#
#  --gtf <string>      finspector.gtf
#
# Optional:
#
#  --kmer_size <int>   default: $KMER_SIZE
#
############################################

__EOUSAGE__

    ;


my $help_flag;
my $fasta_file;
my $gtf_file;


&GetOptions ( 'h' => \$help_flag,
              'fasta=s' => \$fasta_file,
              'gtf=s' => \$gtf_file,
              'kmer_size=i' => \$KMER_SIZE,
    );


if ($help_flag) {
    die $usage;
}

unless ($fasta_file && $gtf_file) {
    die $usage;
}

main: {
    
    my $fasta_reader = new Fasta_reader($fasta_file);
    my %seqs_hash = $fasta_reader->retrieve_all_seqs_hash();
    
    my %contig_to_gene_structs = &parse_gene_pair_info($gtf_file);
    

    print join("\t", "contig", "feature_type", "lend", "rend") . "\n"; # header
    
    foreach my $contig (keys %contig_to_gene_structs) {
        my $contig_seq = uc $seqs_hash{$contig} or die "Error, cannot find seq for $contig";
        
        my ($geneA_struct, $geneB_struct) = @{$contig_to_gene_structs{$contig}};

        # report gene structure info
        my %seen;
        ## include exon structures for vis
        for my $exon_coordset (@{$geneA_struct->{exons}}) {
            my ($exon_lend, $exon_rend) = @$exon_coordset;
            my $token = join("::", "GeneA", $exon_lend, $exon_rend);
            if (! $seen{$token}) {
                print join("\t", $contig, "GeneA", $exon_lend, $exon_rend) . "\n";
                $seen{$token} = 1;
            }
            
        }
        
        for my $exon_coordset (@{$geneB_struct->{exons}}) {
            my ($exon_lend, $exon_rend) = @$exon_coordset;
            my $token = join("::", "GeneA", $exon_lend, $exon_rend);
            if (! $seen{$token}) {
                print join("\t", $contig, "GeneB", $exon_lend, $exon_rend) . "\n";
                $seen{$token} = 1;
            }
        }
                
        my @microhomologies = &find_microhomologies($contig_seq, $geneA_struct, $geneB_struct, $KMER_SIZE);

        if (@microhomologies) {
                    
            foreach my $microhomology (@microhomologies) {
                my ($x, $y) = @$microhomology;
                print join("\t", $contig, "MicroH", $x, $y) . "\n";
            }

        }
    }
    
    exit(0);
}


####
sub parse_gene_pair_info {
    my ($gtf_file) = @_;

    my %gene_to_exon_coords;
    
    open(my $fh, $gtf_file) or die "Error, cannot open file: $gtf_file";
    while(<$fh>) {
        chomp;
        if (/^\#/) { next; }
        unless (/\w/) { next; }
        
        my @x = split(/\t/);

        if ($x[2] eq "exon") {
            my $contig = $x[0];
            my $info = $x[8];
            my $lend = $x[3];
            my $rend = $x[4];
            
            if ($info =~ /gene_id \"([^\"]+)\";/) {
                my $gene_id = $1;
                push (@{$gene_to_exon_coords{$contig}->{$gene_id}}, [$lend, $rend]);
            }
            else {
                die "Error, cannot parse gene info from $info";
            }
        }
    }
    close $fh;

    my %contig_to_gene_structs;

    foreach my $contig (keys %gene_to_exon_coords) {
        
        my $contig_info_href = $gene_to_exon_coords{$contig};
    
        my @gene_ids = keys %$contig_info_href;
        if (scalar(@gene_ids) != 2) {
            die "Error, didn't extract exactly two genes from gtf file";
        }

        my @gene_structs;
        
        foreach my $gene_id (@gene_ids) {
            my @coordsets = @{$contig_info_href->{$gene_id}};
            @coordsets = sort {$a->[0]<=>$b->[0]} @coordsets;
            
            my @exon_coords = @coordsets; # retain original ones before collapsing
            
            @coordsets = &Overlap_piler::simple_coordsets_collapser(@coordsets);
            
            my $range_lend = $coordsets[0]->[0];
            my $range_rend = $coordsets[$#coordsets]->[1];
            
            
            push (@gene_structs, { gene_id => $gene_id,
                                   gene_lend => $range_lend,
                                   gene_rend => $range_rend,
                                   coordsets => \@coordsets,
                                   exons => \@exon_coords,
                  } );
        }
        
        @gene_structs = sort {$a->{gene_lend}<=>$b->{gene_lend}} @gene_structs;

        $contig_to_gene_structs{$contig} = [@gene_structs];
    }
    
    return(%contig_to_gene_structs);
    
}

####
sub find_microhomologies {
    my ($contig_seq, $geneA_struct, $geneB_struct, $KMER_SIZE) = @_;


    ## index kmer positions in geneA
    my %kmer_to_pos_geneA;

    my @exon_coords_geneA = @{$geneA_struct->{coordsets}};

    foreach my $exon_coordset (@exon_coords_geneA) {

        my ($exon_lend, $exon_rend) = @$exon_coordset;

        for (my $i = $exon_lend - $KMER_SIZE+1; $i <= $exon_rend; $i++) {
            my $kmer = substr($contig_seq, $i, $KMER_SIZE);
            push (@{$kmer_to_pos_geneA{$kmer}}, $i);
        }
    }
    
    ## find matches in geneB
    my @exon_coords_geneB = @{$geneB_struct->{coordsets}};

    my @position_hits;
    
    foreach my $exon_coordset (@exon_coords_geneB) {
        
        my ($exon_lend, $exon_rend) = @$exon_coordset;
        
        for (my $i = $exon_lend - $KMER_SIZE+1; $i <= $exon_rend; $i++) {
            my $kmer = substr($contig_seq, $i, $KMER_SIZE);
            
            if (my $posA_hits_aref = $kmer_to_pos_geneA{$kmer}) {
                foreach my $posA_hit (@$posA_hits_aref) {
                    push (@position_hits, [$posA_hit, $i]);
                }
            }
        }
        
    }

    return(@position_hits);
    
}
