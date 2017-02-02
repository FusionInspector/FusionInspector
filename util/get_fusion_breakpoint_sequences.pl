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
use DelimParser;
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);
use Data::Dumper;


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


my $FI_result_file;
my $FI_gtf_file;
my $FI_genome_fasta;

&GetOptions( 'FI_result_file=s' => \$FI_result_file,
             'FI_gtf=s' => \$FI_gtf_file,
             'FI_genome_fasta=s' => \$FI_genome_fasta,
    );

unless ($FI_result_file && $FI_gtf_file && $FI_genome_fasta) {
    die $usage;
}


main: {
    
    my $gene_obj_indexer = {};

    my %fusion_contig_to_coordset_structs;
    
    print STDERR "-parsing GTF file: $FI_gtf_file\n";
    my $asmbl_id_to_gene_list_href = &GTF_utils::index_GTF_gene_objs_from_GTF($FI_gtf_file, $gene_obj_indexer);

    my $fasta_reader = new Fasta_reader($FI_genome_fasta);
    my %genome_seqs = $fasta_reader->retrieve_all_seqs_hash();

    open(my $fh, $FI_result_file) or die "Error, cannot open file $FI_result_file";
    my $tab_reader = new DelimParser::Reader($fh, "\t");

    while (my $row = $tab_reader->get_row()) {
        
        my $fusion_name = $row->{'#FusionName'};

        my $left_breakpoint = $row->{LeftBreakpoint};
        my $right_breakpoint = $row->{RightBreakpoint};

        my $left_local_break = $row->{LeftLocalBreakpoint};
        my $right_local_break = $row->{RightLocalBreakpoint};

        my $contig_seq = $genome_seqs{$fusion_name} or die "Error, could not retrieve contig seq for fusion [$fusion_name]";

        ## create gene graph
        my $coordset_structs_href = $fusion_contig_to_coordset_structs{$fusion_name};
        unless (ref $coordset_structs_href) {
            my @gene_ids = @{$asmbl_id_to_gene_list_href->{$fusion_name}};
            
            $coordset_structs_href = &create_gene_graph_obj($gene_obj_indexer, \@gene_ids);
            
            $fusion_contig_to_coordset_structs{$fusion_name} = $coordset_structs_href;
        }
        
        my @left_coords = &select_exons($coordset_structs_href, $left_local_break, 'left');
        
        my @right_coords = &select_exons($coordset_structs_href, $right_local_break, 'right');
                
        
        unless (@left_coords && @right_coords) {
            print STDERR "-Warning, couldn't identify reference exons overlapping breakpoint for $fusion_name .... skipping.\n";
            next;
        }

        my $left_cdna_seq = &reconstruct_cdna($contig_seq, \@left_coords);
        my $right_cdna_seq = &reconstruct_cdna($contig_seq, \@right_coords);

        print ">$fusion_name|$left_breakpoint|$right_breakpoint\n" . lc($left_cdna_seq) . uc($right_cdna_seq) . "\n";
                
    }
    
    
    exit(0);
}


####
sub reconstruct_cdna {
    my ($contig_seq, $coords_aref) = @_;

    my $seq = "";
    foreach my $coordset (@$coords_aref) {
        my ($lend, $rend) = @$coordset;
        my $region_seq = substr($contig_seq, $lend - 1, $rend - $lend + 1);
        $seq .= $region_seq;
    }

    return($seq);
}


####
sub select_exons {
    my ($coordset_structs_href, $breakpoint, $direction) = @_;

    ## prioritize selected exons by:
    # - exact boundary
    # - overlapping
    # - number isoforms w/ exon
    # - exon length

    my @candidates_exact_boundary;
    my @candidates_overlapping;
    foreach my $struct (values %$coordset_structs_href) {
        
        if ( 
            ($direction eq 'left' && $struct->{rend} == $breakpoint)
            ||
            ($direction eq 'right' && $struct->{lend} == $breakpoint) ) {
            
            push (@candidates_exact_boundary, $struct);
        }
        elsif ($struct->{lend} < $breakpoint && $breakpoint < $struct->{rend}) {
            push (@candidates_overlapping, $struct);
        }
    }
    
    my @final_candidates = (@candidates_exact_boundary) ? @candidates_exact_boundary : @candidates_overlapping;

    unless (@final_candidates) {
        return ();
    }
        
    @final_candidates = sort { $b->{num_isoforms} <=> $a->{num_isoforms}
                               ||
                               $b->{length} <=> $a->{length} } @final_candidates;
    
    my $best_candidate = shift @final_candidates;
    
    my @ret_exons = ($best_candidate);
    my @ret_coords;
    if ($direction eq 'left') {
        
        @ret_coords = ([$best_candidate->{lend}, $breakpoint]);
        while (my $prev_href = $ret_exons[0]->{prev}) {
            my @prev_coords = keys %$prev_href; 
            if (scalar (@prev_coords) == 1) {
                my $struct = $coordset_structs_href->{ $prev_coords[0] };
                unshift(@ret_exons, $struct);
                unshift (@ret_coords, [$struct->{lend}, $struct->{rend}]);
            }
            else {
                last;
            }
        }
    }
    else {
        # right
        @ret_coords = [$breakpoint, $best_candidate->{rend}];
        while (my $next_href = $ret_exons[$#ret_exons]->{next}) {
            my @next_coords = keys %$next_href;
            if (scalar @next_coords == 1) {
                my $struct = $coordset_structs_href->{ $next_coords[0] };
                push (@ret_exons, $struct);
                push (@ret_coords, [$struct->{lend}, $struct->{rend}]);
            }
            else {
                last;
            }
        }
    }

    return(@ret_coords);
}



####
sub create_gene_graph_obj {
    my ($gene_obj_indexer, $gene_ids_aref) = @_;

    my %coordset_to_structs;
    
    foreach my $gene_id (@$gene_ids_aref) {
        my $gene_obj = $gene_obj_indexer->{$gene_id};

        foreach my $transcript_obj ($gene_obj, $gene_obj->get_additional_isoforms()) {
            #print $transcript_obj->toString();
        
            my @exons = $transcript_obj->get_exons();
            
            my $prev_struct = undef;
            foreach my $exon (@exons) {
                my ($lend, $rend) = sort {$a<=>$b} $exon->get_coords();
                
                my $coordstr = "$lend-$rend";
                
                my $struct = $coordset_to_structs{$coordstr};
                unless ($struct) {
                    $struct = $coordset_to_structs{$coordstr} = { 

                        coordstr => $coordstr,
                        
                        lend => $lend,
                        rend => $rend,
                        
                        length => $rend - $lend + 1,
                        
                        prev => {},
                        next => {},
                        
                        num_isoforms => 0,
                    };
                    
                }
                
                $struct->{num_isoforms}++;
                
                if ($prev_struct) {
                    $prev_struct->{next}->{$coordstr}++;
                    $struct->{prev}->{ $prev_struct->{coordstr} }++;
                }
                
                $prev_struct = $struct;
            }
        }
     }

    return(\%coordset_to_structs);

}
