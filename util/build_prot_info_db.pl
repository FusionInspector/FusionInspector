#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);
use Data::Dumper;
use JSON::XS;
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Fasta_reader;
use TiedHash;
use Process_cmd;

my $usage = <<__EOUSAGE__;

##########################################################################
#
#  --ref_coding_GFF3 <string>        : reference annotation for coding genes in gff3 format
#                                      
#  --cds_fasta <string>              : CDS fasta file
#
#  --pfam_domains <string>           : pfam output
#
#  --out_db_dir <string>             : directory to store prot_info.db.idx
#
##########################################################################


 
__EOUSAGE__

    ;


my $help_flag;
my $ref_coding_gff3_file;
my $cds_fasta_file;
my $pfam_domains_file;
my $out_db_dir;

&GetOptions ( 'h' => \$help_flag,
              'ref_coding_GFF3=s' => \$ref_coding_gff3_file,
              'cds_fasta=s' => \$cds_fasta_file,
              'pfam_domains=s' => \$pfam_domains_file,
              'out_db_dir=s' => \$out_db_dir,
    );

unless ($ref_coding_gff3_file && $cds_fasta_file && $pfam_domains_file && $out_db_dir) {
    die $usage;
}



main: {

    unless (-d $out_db_dir) {
        &process_cmd("mkdir -p $out_db_dir");
    }
    my $prot_info_db_idx = "$out_db_dir/prot_info_db.idx";
    if (-e $prot_info_db_idx) {
        die "Error, $prot_info_db_idx already exists.  Please remove it or rename it before proceeding";
    }
        
    my $annot_manager = Annotation_manager->new($ref_coding_gff3_file);

    print STDERR "-parsing cds fa: $cds_fasta_file\n";
    my $fasta_reader = new Fasta_reader($cds_fasta_file);
    my %cds_seqs = $fasta_reader->retrieve_all_seqs_hash();
    
    my %pfam_hits = &parse_pfam($pfam_domains_file);
    
    $annot_manager->add_cds_and_pfam(\%cds_seqs, \%pfam_hits);
    
    $annot_manager->build_prot_info_db($prot_info_db_idx);
    
    
    exit(0);
}


=pfam_entry

0       BTG
1       PF07742.9
2       116
3       BTG2|ENST00000290551.4
4       -
5       158
6       3.8e-42
7       142.7
8       0.0
9       1
10      1
11      2.8e-46
12      4.5e-42
13      142.4
14      0.0
15      1
16      115
17      9
18      121
19      9
20      122
21      0.99
22      BTG
23      family

=cut


####
sub parse_pfam {
    my ($pfam_file) = @_;

    print STDERR "-parsing $pfam_file\n";
    
    my %model_to_domains;
    
    open (my $fh, $pfam_file) or die "Error, cannot open file $pfam_file";
    while (<$fh>) {
        chomp;
        unless (/\w/) { next; }
        if (/^\#/) { next; }
        my @x = split(/\s+/);
        
        if (scalar @x < 22) {
            print STDERR "WARNING: Skipping line: $_ as likely corrupt.\n";
            next;
        }
        
        my $QueryProtID = $x[3];
        my $pfam_id = $x[1];
        my $HMMERDomain = $x[0];
        my $HMMERTDomainDescription = join(" ", @x[22..$#x]);
        my $QueryStartAlign = $x[17];
        my $QueryEndAlign = $x[18];
        my $PFAMStartAlign = $x[15];
        my $PFAMEndAlign = $x[16];
        my $FullSeqEvalue = $x[6];
        my $ThisDomainEvalue = $x[11];
        my $FullSeqScore = $x[7];
        my $FullDomainScore = $x[13];
        
        my $pfam_domain_struct = { cds_id => $QueryProtID,
                                   hmmer_domain => $HMMERDomain,
                                   query_start => $QueryStartAlign,
                                   query_end => $QueryEndAlign,
                                   domain_evalue => $ThisDomainEvalue,
        };
        
        push (@{$model_to_domains{$QueryProtID}}, $pfam_domain_struct);
    }
    close $fh;

    return(%model_to_domains);
}
        
        

#############################################################################
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


sub TO_JSON {
    return { %{ shift() } };
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
sub add_cds_and_pfam {
    my ($self) = shift;
    my ($cds_seqs_href, $pfam_hits_href) = @_;
    
    
    my @gene_ids = $self->get_gene_list();
    
    foreach my $gene_id (@gene_ids) {
        my @cds_features = $self->get_CDS_features($gene_id);
        foreach my $cds_feature (@cds_features) {

            my $cds_id = $cds_feature->{cds_id};
            my $cds_seq = $cds_seqs_href->{$cds_id};
            unless ($cds_seq) {
                print STDERR "WARNING, no CDS sequence for $cds_id\n";
                $cds_seq = "";
            }
            
            $cds_feature->set_CDS_sequence($cds_seq);

            my $pfam_hits = $pfam_hits_href->{$cds_id};
            if (ref $pfam_hits) {
                $cds_feature->add_pfam_hits(@$pfam_hits);
            }
        }
    }

    return;
}

####
sub build_prot_info_db {
    my ($self) = shift;
    my ($prot_db_idx_file) = @_;
    
    my $tied_hash = new TiedHash( { create => $prot_db_idx_file } );
    
    
    my @gene_ids = $self->get_gene_list();
    
    my $coder = JSON::XS->new->convert_blessed;
    
    foreach my $gene_id (@gene_ids) {
        my @cds_features = $self->get_CDS_features($gene_id);
        

        my $gene_id_store = $gene_id;
        $gene_id_store =~ s/\|.*$//;

        my $json = $coder->pretty->encode(\@cds_features);

        $tied_hash->store_key_value($gene_id_store, $json);
    }
    
    $tied_hash = undef; # closes it.
    
    return;
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
        my $chr = $x[0];
        
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
            
            $self->_add_to_CDS_feature($gene_id, $mRNA_id, $chr, $start, $end, $orient, $phase_init);
        }
    }

    $self->_refine_CDS_features();
}

####
sub _add_to_CDS_feature {
    my ($self) = shift;
    my ($gene_id, $cds_id, $chr, $start, $end, $orient, $phase_init) = @_;

    my $cds_feature_obj = $self->_get_CDS_feature($gene_id, $cds_id);
    
    $cds_feature_obj->add_segment($chr, $start, $end, $orient, $phase_init);
    
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
        
        refined_flag => 0,

        pfam_hits => [],

        cds_seq => "",
        
        
    };

    bless($self, $packagename);

    return($self);
}


sub TO_JSON {
    return { %{ shift() } };
}


####
sub add_segment {
    my ($self) = shift;
    my ($chr, $lend, $rend, $orient, $phase_beg) = @_;

    my $phased_segment = { chr => $chr,
                           
                           lend => $lend,
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

        my $phase_end = ".";
        if ($phase_beg ne ".") {
            my $adj_seg_len = $seg_len;
            $adj_seg_len += $phase_beg;
        
            $phase_end = ($adj_seg_len -1)  % 3;
        }
        
        $segment->{rel_lend} = $rel_lend;
        $segment->{rel_rend} = $rel_rend;
        $segment->{phase_end} = $phase_end;
        
        $sum_segs_len += $seg_len;
    
    }

    $self->{refined_flag} = 1;

    return;
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

####
sub set_CDS_sequence {
    my ($self) = shift;
    my ($cds_seq) = @_;

    $self->{cds_seq} = $cds_seq;

    return;
}

####
sub add_pfam_hits {
    my $self = shift;
    my (@pfam_hits) = @_;

    push (@{$self->{pfam_hits}}, @pfam_hits);

    return;
}
