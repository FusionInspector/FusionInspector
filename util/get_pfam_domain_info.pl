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
##########################################################################


__EOUSAGE__

    ;


my $help_flag;
my $finspector_gtf;
my $genome_lib_dir;


&GetOptions ( 'h' => \$help_flag,
              'finspector_gtf=s' => \$finspector_gtf,
              'genome_lib_dir=s' => \$genome_lib_dir,
    );

unless ($finspector_gtf && $genome_lib_dir) {
    
    die $usage;
}


main: {

    my $pfam_domain_dbm = "$genome_lib_dir/pfam_domains.dbm";
    if (! -s $pfam_domain_dbm) {
        die "Error, cannot locate pfam ctat database: $pfam_domain_dbm";
    }
    
    my $pfam_domain_db_tied_hash = new TiedHash( { 'use' => $pfam_domain_dbm } );
    
    my %cds_to_finspector_coordinates = &parse_CDS_info($finspector_gtf);
    

    my %contig_to_pfam_coords;

    foreach my $contig (keys %cds_to_finspector_coordinates) {

        my $contig_trans_href = $cds_to_finspector_coordinates{$contig};
        
        my @transcripts = keys %{$contig_trans_href};
        
        foreach my $transcript (@transcripts) {
                        
            if (my @pfam_domains = &get_pfam_domains($transcript, $pfam_domain_db_tied_hash)) {
                
                ## adjust pfam domains for it

                my @cds_coords_info = @{$contig_trans_href->{$transcript}};
                
                &assign_rel_coords(@cds_coords_info);
                
                foreach my $pfam_domain_struct (@pfam_domains) {
                    
                    &capture_contig_pfam_domain(\%contig_to_pfam_coords, $contig, \@cds_coords_info, $pfam_domain_struct);
                    
                }
                

            }
        }
        
                    
    }
    
    ## Resolve overlaps, write output coordinates.
    foreach my $contig (keys %contig_to_pfam_coords) {
        my $pfam_entries_href = $contig_to_pfam_coords{$contig};
        
        my @pfam_accs = keys %$pfam_entries_href;
        foreach my $pfam_acc (@pfam_accs) {
            my @coordsets = @{$pfam_entries_href->{$pfam_acc}};
            my @collapsed_coordsets = &Overlap_piler::simple_coordsets_collapser(@coordsets);
            
            foreach my $collapsed_coordset (@collapsed_coordsets) {
                my ($match_lend, $match_rend) = @$collapsed_coordset;

                print join("\t", $contig, "Pfam", "match", $match_lend, $match_rend, ".", "+", ".", "$pfam_acc") . "\n";
            }
        }
    }
    
    exit(0);
    
}

####
sub parse_CDS_info {
    my ($finspector_gtf) = @_;

    my %cds_to_finspector_coordinates;

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
        
        if ($feat_type eq "CDS") {
            # transcript_id "AC099850.1--VMP1^ENST00000588915.1"
            if ($info =~ /transcript_id \"([^\"]+)\"/) {
                my $transcript_id = $1;
                my @pts = split(/\^/, $transcript_id);
                $transcript_id = pop @pts;
            
                push (@{$cds_to_finspector_coordinates{$fusion_contig}->{$transcript_id}}, 
                      { 
                          contig_coords => [ $lend, $rend ] 
                      } );
            }
        }
    }
    close $fh;
    
    return(%cds_to_finspector_coordinates);
    
}



####
sub get_pfam_domains {
    my ($cds_id, $pfam_domain_db_tied_hash) = @_;

    unless ($pfam_domain_db_tied_hash) {
        return;
    }

    my $pfam_hits_aref = &decode_my_json($cds_id, $pfam_domain_db_tied_hash);
    
    if ($pfam_hits_aref) {
        return(@$pfam_hits_aref);
    }
    else {
        return ();
    }
                
}



sub decode_my_json {
    my ($key, $prot_info_db_tied_hash, $report_failed_retrievals_flag) = @_;
    
    my $coder = new JSON::XS();

    my $decoded = undef;
    my $json = undef;
    
    eval {
        $json = $prot_info_db_tied_hash->get_value($key);

        if ($json) {
            $decoded = $coder->decode($json);
            #print Dumper($decoded);
        }
        else {
            if ($report_failed_retrievals_flag) {
                print STDERR "WARNING, no entry stored in dbm for [$key]\n";
            }
            return(undef);
        }
    };

    if ($@) {
        print STDERR "WARNING, key: $key returns json: $json and error decoding: $@";
    }

    return($decoded);
}


####
sub assign_rel_coords {
    my @cds_coords_info = @_;

    @cds_coords_info = sort {$a->{contig_coords}->[0] <=> $b->{contig_coords}->[0] } @cds_coords_info;

    my $prev_lend = 0;
    foreach my $cds_coords_struct (@cds_coords_info) {
        my ($contig_lend, $contig_rend) = @{$cds_coords_struct->{contig_coords}};
        
        my $rel_lend = $prev_lend + 1;
        my $rel_rend = $prev_lend + ($contig_rend - $contig_lend + 1);
        
        $prev_lend = $rel_rend;
        $cds_coords_struct->{rel_coords} = [$rel_lend, $rel_rend];
    }

    return;
}


####
sub capture_contig_pfam_domain {
    my ($contig_to_pfam_coords_href, $contig, $cds_coords_info_aref, $pfam_domain_struct) = @_;

    my $pfam_domain = $pfam_domain_struct->{hmmer_domain};
    my $pfam_lend = $pfam_domain_struct->{query_start} * 3; # nts instead of protein coords
    my $pfam_rend = $pfam_domain_struct->{query_end} * 3;

    foreach my $coordset (@$cds_coords_info_aref) {
        
        my ($contig_lend, $contig_rend) = @{$coordset->{contig_coords}};
        my ($rel_lend, $rel_rend) = @{$coordset->{rel_coords}};

        if ($pfam_lend < $rel_rend && $pfam_rend > $rel_lend) {
            # overlap!
            my $pfam_end5 = ($pfam_lend <= $rel_lend) ? $contig_lend : $pfam_lend - $rel_lend + $contig_lend;
            my $pfam_end3 = ($pfam_rend >= $rel_rend) ? $contig_rend : $contig_rend - ($rel_rend - $pfam_rend);

            
            push (@{$contig_to_pfam_coords_href->{$contig}->{$pfam_domain}}, [$pfam_end5, $pfam_end3] );
        }
    }

    return;


}
