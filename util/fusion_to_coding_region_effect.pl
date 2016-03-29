#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);
use Data::Dumper;


my $usage = <<__EOUSAGE__;

##########################################################################
#
#  --fusions <string>                : fusion predictions
#
#  --ref_coding_GFF3 <string>        : reference annotation for coding genes in gff3 format
#                                      
#  --pfam_domains <string>           : pfam output
#
##########################################################################


__EOUSAGE__

    ;


my $help_flag;
my $fusions_file;
my $ref_coding_gff3_file;
my $pfam_domains_file;

&GetOptions ( 'h' => \$help_flag,
              'fusions=s' => \$fusions_file,
              'ref_coding_GFF3=s' => \$ref_coding_gff3_file,
              'pfam_domains=s' => \$pfam_domains_file,
              
    );

unless ($fusions_file && $ref_coding_gff3_file && $pfam_domains_file) {
    die $usage;
}


main: {

    #my @fusions = &parse_fusions($fusions_file);

    my $annot_manager = Annotation_manager->new($ref_coding_gff3_file);

    #print Dumper($annot_manager);

    print $annot_manager->toString();
    

    exit(0);
}

####
package Annotation_manager;
use strict;
use warnings;

####
sub new {
    my ($packagename, $gff3_file) = @_;

    my $self = {
        gene_to_CDS_features => {},  # gene_id -> cds_id -> cds_feature_obj
    };
    
    bless($self, $packagename);

    $self->parse_GFF3_instantiate_featureset($gff3_file);

    
    
    return($self);
}


####
sub get_gene_list {
    my ($self) = @_;

    my @gene_ids = keys %{$self->{gene_to_CDS_features}};

    return(@gene_ids);
}

####
sub get_CDS_features {
    my $self = shift;
    my ($gene_id) = @_;

    my $cds_features_href = $self->{gene_to_CDS_features}->{$gene_id};
    
    if (ref $cds_features_href) {
        return(values %$cds_features_href);
    }

    else {
        return();
    }
}

####
sub toString {
    my ($self) = shift;
    
    my @gene_ids = $self->get_gene_list();

    my $text = "";
    
    foreach my $gene_id (@gene_ids) {
        my @cds_features = $self->get_CDS_features($gene_id);
        foreach my $cds_feature (@cds_features) {
            $text .= $cds_feature->toString() . "\n";
        }
    }
    
    return ($text);
}
    

    
####
sub parse_GFF3_instantiate_featureset {
    my ($self) = shift;
    my ($gff3_file) = @_;
    
    my %mRNA_to_gene;

    print STDERR "-parsing gff3 coding regions file: $gff3_file\n";
    open(my $fh, $gff3_file) or die "Error, cannot open file $gff3_file";
    while (<$fh>) {
        chomp;
        unless (/\w/) { next; }
        if (/^\#/) { next; }
        
        my @x = split(/\t/);

        my %info = &_parse_info($x[8]);

        if ($x[2] eq 'mRNA') {
            my $mRNA_id = $info{ID};
            my $gene_id = $info{Parent};
            $mRNA_to_gene{$mRNA_id} = $gene_id;
        }
        elsif ($x[2] eq 'CDS') {
            my $cds_id = $info{ID};
            my $mRNA_id = $info{Parent};
            my $gene_id = $mRNA_to_gene{$mRNA_id} or die "Error, no gene_id for mRNA: $mRNA_id";
            
            my $start = $x[3];
            my $end = $x[4];
            my $orient = $x[6];
            my $phase_init = $x[7];
            
            $phase_init = &_phase_transform($phase_init);
            
            $self->_add_to_CDS_feature($gene_id, $mRNA_id, $start, $end, $orient, $phase_init);
        }
    }

    $self->_refine_CDS_features();
}

####
sub _add_to_CDS_feature {
    my ($self) = shift;
    my ($gene_id, $cds_id, $start, $end, $orient, $phase_init) = @_;

    my $cds_feature_obj = $self->_get_CDS_feature($gene_id, $cds_id);
    
    $cds_feature_obj->add_segment($start, $end, $orient, $phase_init);
    
}

####
sub _get_CDS_feature {
    my ($self) = shift;
    my ($gene_id, $cds_id) = @_;

    my $cds_feature_obj = $self->{gene_to_CDS_features}->{$gene_id}->{$cds_id};

    unless (ref $cds_feature_obj) {
        $cds_feature_obj = $self->{gene_to_CDS_features}->{$gene_id}->{$cds_id} = CDS_feature->new($gene_id, $cds_id);
    }

    return($cds_feature_obj);
}


####
sub _phase_transform {
    my ($phase) = @_;

    if ($phase eq "0") {
        return(0);
    }
    elsif ($phase eq "1") {
        return(2);
    }
    elsif ($phase eq "2") {
        return(1);
    }
    else {
        return($phase);
    }
}

    
####
sub _parse_info {
    my ($info) = @_;

    my %ret;

    my @keyvals = split(/;/, $info);
    foreach my $keyval (@keyvals) {
        my ($key, $val) = split(/=/, $keyval);
        $ret{$key} = $val;
    }
    
    return(%ret);
}


####
sub _refine_CDS_features {
    my $self = shift;

    my @gene_ids = $self->get_gene_list();

    foreach my $gene_id (@gene_ids) {
        
        my $cds_features_href = $self->{gene_to_CDS_features}->{$gene_id};
        my @cds_feature_objs = values %$cds_features_href;

        foreach my $cds_feature_obj (@cds_feature_objs) {
            $cds_feature_obj->refine();
        }
    }
    
    return;
}



####################################################################################
package CDS_feature;
use strict;
use warnings;

####
sub new {
    my ($packagename) = shift;
    my ($gene_id, $cds_id) = @_;

    
    my $self = {
        phased_segments => [],
        
        gene_id => $gene_id,
        cds_id => $cds_id,
                
        init_flag => 0,
    
    };

    bless($self, $packagename);

    return($self);
}

####
sub add_segment {
    my ($self) = shift;
    my ($lend, $rend, $orient, $phase_beg) = @_;

    my $phased_segment = { lend => $lend,
                           rend => $rend,
                           orient => $orient,
                           phase_beg => $phase_beg,

                           rel_lend => undef,
                           rel_rend => undef,
                           phase_end => undef,  ## all set on init
    };

    push (@{$self->{phased_segments}}, $phased_segment);
    
    return;
}

####
sub get_segments {
    my ($self) = shift;

    return(@{$self->{phased_segments}});

}


####
sub refine {
    my ($self) = shift;
    
    my @segments = $self->get_segments();

    @segments = sort {$a->{lend} <=> $b->{lend}} @segments;
    
    my $orient = $segments[0]->{orient};

    if ($orient eq '-') {
        @segments = reverse @segments;
    }

    my $sum_segs_len = 0;
    foreach my $segment (@segments) {

        my $seg_len = $segment->{rend} - $segment->{lend} + 1;
        my $phase_beg = $segment->{phase_beg};

        my $rel_lend = $sum_segs_len + 1;
        my $rel_rend = $sum_segs_len + $seg_len;
     
        my $adj_seg_len = $seg_len;
        $adj_seg_len += $phase_beg;
        
        my $phase_end = ($adj_seg_len -1)  % 3;
        
        $segment->{rel_lend} = $rel_lend;
        $segment->{rel_rend} = $rel_rend;
        $segment->{phase_end} = $phase_end;
        
        $sum_segs_len += $seg_len;
    
    }
        
}


####
sub toString {
    my ($self) = shift;
    
    my @segments = $self->get_segments();
    
    @segments = sort {$a->{lend} <=> $b->{lend}} @segments;

    my $orient = $segments[0]->{orient};
    if ($orient eq '-') {
        @segments = reverse @segments;
    }


    my $ret_text = "";
    
    foreach my $segment (@segments) {
        
        $ret_text .= join("\t", $self->{gene_id}, $self->{cds_id}, 
                          $segment->{lend}, $segment->{rend},
                          $segment->{orient}, 
                          $segment->{rel_lend}, $segment->{rel_rend},
                          $segment->{phase_beg}, $segment->{phase_end}) . "\n";
    }
    
    
    
    return ($ret_text);
        
}
