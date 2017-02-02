#!/usr/bin/env perl

# borrowed from PASA/misc_utilities

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Gene_obj;
use Gene_obj_indexer;
use GTF_utils;
use Fasta_reader;
use Carp;
use Delim_parser;
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);


my $usage = <<__EOUSAGE__;

#######################################################################
#
#  --FI_result_file <string>       FusionInspector result file
#
#  --FI_gtf <string>               FusionInspector-generated gtf file
#
#  --FI_genome_fasta <string>      FusionInspector-generated fasta file
#
#######################################################################


__EOUSAGE__

    ;


&GetOptions( 'FI_result_file=s' => \$FI_result_file,
             'FI_gtf=s' => \$FI_gtf_file,
             'FI_genome_fasta=s' => \$FI_genome_fasta,
    );


main: {
    
    my $gene_obj_indexer = {};

    my %gene_graph_objs;
    
    print STDERR "-parsing GTF file: $gtf_file\n";
    my $asmbl_id_to_gene_list_href = &GTF_utils::index_GTF_gene_objs_from_GTF($FI_gtf_file, $gene_obj_indexer);

    my $fasta_reader = new Fasta_reader($FI_genome_fasta);
    my %genome_seqs = $fasta_reader->retrieve_all_seqs_hash();

    open(my $fh, $FI_result_file) or die "Error, cannot open file $FI_result_file";
    my $tab_reader = new Delim_parser($fh, "\t");

    while (my $row = $tab_reader->get_row()) {
        
        my $fusion_name = $row->{'#FusionName'};

        my $left_breakpoint = $row->{LeftBreakpoint};
        my $right_breakpoint = $row->{RightBreakpoint};

        my $left_local_break = $row->{LeftLocalBreakpoint};
        my $right_local_break = $row->{RightLocalBreakpoint};

        my $contig_seq = $genome_seqs{$fusion_name} or die "Error, could not retrieve contig seq for fusion [$fusion_name]";

        ## create gene graph
        my $gene_graph_obj = $gene_graph_objs{$fusion_name};
        unless (ref $gene_graph_obj) {
            my @gene_ids = @{$asmbl_id_to_gene_list_href->{$fusion_name}};
            
            $gene_graph_obj = &create_gene_graph_obj($gene_obj_indexer, \@gene_ids);
        }
                
        
        my @left_exons = &select_exon($gene_obj_indexer, $gene_graph, $left_local_breakpoint, 'left');
        my @right_exons = &select_exon($gene_obj_indexer, $gene_graph, $right_local_breakpint, 'right');

        
        
    }
            
    
    exit(0);
}

####
sub create_gene_graph_obj {
    my ($gene_obj_indexer, $gene_ids_aref) = @_;

    my %coordset_to_exons;
    
    foreach my $gene_id (@gene_ids_aref) {
        my $gene_obj = $gene_obj_indexer->{$gene_id};

        foreach my $transcript_obj ($gene_obj, $gene_obj->get_additional_isoforms()) {

        }
    }

}
