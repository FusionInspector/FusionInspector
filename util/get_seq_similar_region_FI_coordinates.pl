#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);
use Data::Dumper;
use JSON::XS;
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use TiedHash;
use Overlap_piler;

my $DEBUG = 0;

my $usage = <<__EOUSAGE__;

##########################################################################
#
#  --finspector_gtf <string>                : fusion predictions
#
#  --genome_lib_dir <string>                : CTAT genome lib
#            
#  --debug|d                                debug mode
#                          
##########################################################################


__EOUSAGE__

    ;


my $help_flag;
my $finspector_gtf;
my $genome_lib_dir;


&GetOptions ( 'h' => \$help_flag,
              'finspector_gtf=s' => \$finspector_gtf,
              'genome_lib_dir=s' => \$genome_lib_dir,
    
              'debug|d' => \$DEBUG,
    );

unless ($finspector_gtf && $genome_lib_dir) {
    
    die $usage;
}


my $BLAST_ALIGNS_IDX;
my $blast_aligns_idx_file = "$genome_lib_dir/trans.blast.align_coords.align_coords.dbm";
if (-s $blast_aligns_idx_file) {
    $BLAST_ALIGNS_IDX = new TiedHash( { use => $blast_aligns_idx_file } );
}
else {
    die "Error, cannot lcoate blast idx file: $blast_aligns_idx_file";
}

my $JSON_DECODER = JSON::XS->new();

my $MATCH_COUNTER;

main: {
    
    my %cdna_to_finspector_coordinates = &parse_cdna_info($finspector_gtf);
    
    foreach my $contig (keys %cdna_to_finspector_coordinates) {
        
        my $blast_pair_info = $BLAST_ALIGNS_IDX->get_value($contig);
        
        unless ($blast_pair_info) { next; }
        
        my $blast_align_info_struct = $JSON_DECODER->decode($blast_pair_info);
    
        my ($geneA, $geneB) = split(/--/, $contig);
        
        my $coord_mappings_href = $cdna_to_finspector_coordinates{$contig};
        print Dumper($coord_mappings_href) if $DEBUG;

        my $seq_similar_coords_A_aref = $blast_align_info_struct->{coords_A};
        my $FI_coord_mappings_A_aref = $coord_mappings_href->{$geneA};

        if ($DEBUG) {
            print "geneA $geneA seq_similar_coords: " . Dumper($seq_similar_coords_A_aref);
            print "FI coord mappings for $geneA: " . Dumper($FI_coord_mappings_A_aref);
        }
        
        my @seq_similar_geneA_regions = &define_seq_similar_FI_regions($seq_similar_coords_A_aref, $FI_coord_mappings_A_aref);
        
        &write_gtf($contig, ++$MATCH_COUNTER, \@seq_similar_geneA_regions);
        
        my $seq_similar_coords_B_aref = $blast_align_info_struct->{coords_B};
        my $FI_coord_mappings_B_aref = $coord_mappings_href->{$geneB};
        
        my @seq_similar_geneB_regions = &define_seq_similar_FI_regions($seq_similar_coords_B_aref, $FI_coord_mappings_B_aref);
                
        
        &write_gtf($contig, ++$MATCH_COUNTER, \@seq_similar_geneB_regions);
        
        
    }
    

    exit(0);

}


####
sub write_gtf {
    my ($contig, $match_count, $seq_similar_regions_aref) = @_;

    foreach my $collapsed_coordset (@$seq_similar_regions_aref) {
        my ($match_lend, $match_rend) = @$collapsed_coordset;
        
        print join("\t", $contig, "seqsimilar", "match", $match_lend, $match_rend, ".", "+", ".", "ID=seqsimilar_$MATCH_COUNTER; Target=blastpair") . "\n";
    }
    
    return;
        
}

####
sub parse_cdna_info {
    my ($finspector_gtf) = @_;

    my %cdna_to_finspector_coordinates;

    my %seen;

    open(my $fh, $finspector_gtf) or die "Error, cannot open file: $finspector_gtf";
    while(<$fh>) {
        chomp;
        unless (/\w/) { next; }
        my @x = split(/\t/);
        
        my $fusion_contig = $x[0];
        my $feat_type = $x[2];
        my $lend = $x[3];
        my $rend = $x[4];
        
        my $info = $x[8];
        
        if ($feat_type eq "exon") {
            # transcript_id "AC099850.1--VMP1^ENST00000588915.1"
            
            my $gene_name;
            my $orig_coords;

            
            if ($info =~ /gene_name \"([^\"]+)\"/) {
                $gene_name = $1;
            }
            elsif ($info =~ /gene_id \"([^\"]+)/) {
                $gene_name = $1;
            }
            else {
                die "Error, couldn't extract gene_name or gene_id from $info";
            }
            
            if ($info =~ /orig_coord_info \"([^\"]+)\"/) {
                $orig_coords = $1;
            }
            else {
                die "Error, couldn't extract orig coords from $info";
            }
            
            my $token = join("^", $fusion_contig, $gene_name, $lend, $rend, $orig_coords);
            unless ($seen{$token}) {

                my ($chr, $genome_lend, $genome_rend, $genome_orient) = split(/,/, $orig_coords);


                push (@{$cdna_to_finspector_coordinates{$fusion_contig}->{$gene_name}}, 
                      { 
                          contig_coords => [ $lend, $rend ],
                          genome_coords => [ $genome_lend, $genome_rend],
                          genome_orient => $genome_orient,,
                      } );
                
                $seen{$token} = 1;
            }
        }
        
    }
    close $fh;
    
        
    return(%cdna_to_finspector_coordinates);
    
}



####
sub assign_rel_coords {
    my @cdna_coords_info = @_;

    @cdna_coords_info = sort {$a->{contig_coords}->[0] <=> $b->{contig_coords}->[0] } @cdna_coords_info;

    my $prev_lend = 0;
    foreach my $cdna_coords_struct (@cdna_coords_info) {
        my ($contig_lend, $contig_rend) = @{$cdna_coords_struct->{contig_coords}};
        
        my $rel_lend = $prev_lend + 1;
        my $rel_rend = $prev_lend + ($contig_rend - $contig_lend + 1);
        
        $prev_lend = $rel_rend;
        $cdna_coords_struct->{rel_coords} = [$rel_lend, $rel_rend];
    }

    return;
}

####
sub define_seq_similar_FI_regions {
    my ($seq_sim_coords_aref, $FI_coord_mappings_aref) = @_;

    my @coordpairs;

    foreach my $seq_sim_coordset (@$seq_sim_coords_aref) {
        my ($seq_sim_genome_lend, $seq_sim_genome_rend) = sort {$a<=>$b} @$seq_sim_coordset;
        
        foreach my $FI_coord_mapping (@$FI_coord_mappings_aref) {
            my ($contig_lend, $contig_rend) = sort {$a<=>$b} @{$FI_coord_mapping->{contig_coords}};

            my ($genome_lend, $genome_rend) = sort {$a<=>$b} @{$FI_coord_mapping->{genome_coords}};

            my $genome_orient = $FI_coord_mapping->{genome_orient};

            
            if ($genome_lend < $seq_sim_genome_rend && $genome_rend > $seq_sim_genome_lend) {
                # overlap!

                my @genome_coords = sort {$a<=>$b} ($genome_lend, $seq_sim_genome_lend, $genome_rend, $seq_sim_genome_rend);
                
                my ($overlap_genome_lend, $overlap_genome_rend) = ($genome_coords[1], $genome_coords[2]);
                

                my $lend_delta = $overlap_genome_lend - $genome_lend;
                my $rend_delta = $genome_rend - $overlap_genome_rend;
                
                if ($genome_orient eq '+') {
                    
                    #                 sim_lend |--| seq_sim_rend
                    #           genome_lend |--------| genome_rend
                    #           (contig_lend)          (contig_rend)
                    
                    my $seq_sim_contig_lend = $contig_lend + $lend_delta;
                    
                    my $seq_sim_contig_rend  = $contig_rend - $rend_delta;

                    push (@coordpairs, [$seq_sim_contig_lend, $seq_sim_contig_rend]);
                }
                else {
                    # '-' genome orient
                    
                    #                 sim_lend |--| sim_rend
                    #           genome_lend |--------| genome_rend
                    #           (contig_rend)          (contig_lend)
                    
                    my $seq_sim_contig_lend = $contig_lend + $rend_delta;
                    
                    my $seq_sim_contig_rend  = $contig_rend - $lend_delta;

                    push (@coordpairs, [$seq_sim_contig_lend, $seq_sim_contig_rend]);
                }
            }
        }

    }


    ## collapse seq-similar regions.
    my @collapsed_coords = &Overlap_piler::simple_coordsets_collapser(@coordpairs);
    
    return(@collapsed_coords);
}
