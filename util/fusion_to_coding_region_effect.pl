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
use Nuc_translator;

my $usage = <<__EOUSAGE__;

##########################################################################
#
#  --fusions <string>                : fusion predictions
#
#  --prot_info_db <string>           :  prot_info_db.idx - path to file.
#                                      
##########################################################################


__EOUSAGE__

    ;


my $help_flag;
my $fusions_file;
my $prot_info_db;

&GetOptions ( 'h' => \$help_flag,
              'fusions=s' => \$fusions_file,
              'prot_info_db=s' => \$prot_info_db,
    );

unless ($fusions_file && $prot_info_db) {
    
    die $usage;
}


main: {

    open (my $fh, $fusions_file) or die "Error, cannot open file $fusions_file";
    my $header = <$fh>;
    unless ($header =~ /^\#fusion_name/) {
        die "Error, fusion_file: $fusions_file has unrecognizable header: $header";
    }
    print $header;

    
    my $tied_hash = new TiedHash( { 'use' => $prot_info_db } );

    while (<$fh>) {
        
        my $line = $_;
        chomp;
        my @x = split(/\t/);

        my $fusion_name = $x[0];
        my $gene_left = $x[6];
        my $break_left = $x[7];
                
        my $gene_right = $x[8];
        my $break_right = $x[9];

        my @results = &examine_fusion_coding_effect($gene_left, $break_left, $gene_right, $break_right, $tied_hash);
        
        foreach my $result (@results) {
            print join("\t", $fusion_name, 
                       $result->{cds_left_id}, $result->{cds_left_range},
                       $result->{cds_right_id}, $result->{cds_right_range},
                       $result->{prot_fusion_type},
                       $result->{fusion_coding_descr},
                       $result->{cds_fusion_seq},
                       $result->{prot_fusion_seq} ) . "\n";
        }
    }
    close $fh;
    
    #my $coder = JSON::XS->new->convert_blessed;
    #print $coder->pretty->encode($annot_manager);
    
    
    exit(0);
}


####
sub examine_fusion_coding_effect {
    my ($gene_left, $break_left, $gene_right, $break_right, $tied_hash) = @_;

    my $coder = new JSON::XS();
        
    my $gene_left_json = $tied_hash->get_value($gene_left);    
    my $gene_left_aref = $coder->decode($gene_left_json);

    print Dumper($gene_left_aref);
    #print "Gene_left: $gene_left_json\n";
    
    
    my $gene_right_json = $tied_hash->get_value($gene_right);
    my $gene_right_aref = $coder->decode($gene_right_json);

    #print Dumper($gene_right_aref);
    #print "Gene_right: $gene_right_json\n";

    my @results;

    
    foreach my $cds_left_obj (@$gene_left_aref) {
        
        my $cds_left_seq = $cds_left_obj->{cds_seq};
        my $cds_left_id = $cds_left_obj->{cds_id};
        
        foreach my $cds_right_obj (@$gene_right_aref) {
            
            my $cds_right_seq = $cds_right_obj->{cds_seq};
            my $cds_right_id = $cds_right_obj->{cds_id};
            
            my ($left_fuse_segments_aref, $right_fuse_segments_aref) = &try_fuse_cdss($cds_left_obj, $break_left, $cds_right_obj, $break_right);

            if (@$left_fuse_segments_aref && @$right_fuse_segments_aref) {

                ## see if compatible
                my $terminal_left_seg = $left_fuse_segments_aref->[$#$left_fuse_segments_aref];
                my $left_end_phase = $terminal_left_seg->{phase_end};
                my $left_rel_rend = $terminal_left_seg->{rel_rend};

                my $left_cds_part = substr($cds_left_seq, 0, $left_rel_rend);
                
                my $initial_right_seg = $right_fuse_segments_aref->[0];
                my $right_beg_phase = $initial_right_seg->{phase_beg};
                my $right_rel_lend = $initial_right_seg->{rel_lend};
                
                my $right_cds_part = substr($cds_right_seq, $right_rel_lend - 1);

                my $fusion_seq = join("", lc($left_cds_part), uc($right_cds_part));
                my $pep = translate_sequence($fusion_seq, 1);

                my $prot_fusion_type = ( ($left_end_phase + 1) % 3 == $right_beg_phase) ? "INFRAME" : "FRAMESHIFT";
                
                my $left_segs_string = &segments_to_string(@$left_fuse_segments_aref);
                my $right_segs_string = &segments_to_string(@$right_fuse_segments_aref);
                
                push (@results, { cds_left_id => $cds_left_id,
                                  cds_right_id => $cds_right_id,
                                  cds_left_range => "1-$left_rel_rend",
                                  cds_right_range => "$right_rel_lend-" . length($cds_right_seq),
                                  prot_fusion_type => $prot_fusion_type,
                                  cds_fusion_seq => $fusion_seq,
                                  prot_fusion_seq => $pep,
                                  fusion_coding_descr => join("<==>", $left_segs_string, $right_segs_string),
                      }
                    );
                
            }
        }
    }
    
    return (@results);
}


####
sub try_fuse_cdss {
    my ($cds_left_obj, $break_left, $cds_right_obj, $break_right) = @_;

    
    # get left part
    my @left_fusion_partner_segments = &get_left_fusion_partner_segments($cds_left_obj, $break_left);
        
    # todo: get right part
    my @right_fusion_partner_segments = &get_right_fusion_partner_segments($cds_right_obj, $break_right);
    
        
    ## piece it together.
    
    #print STDERR "Left: " . Dumper(\@left_fusion_partner_segments) . "\nRight: " . Dumper(\@right_fusion_partner_segments);
    
    return (\@left_fusion_partner_segments, \@right_fusion_partner_segments);
    
}


####
sub get_left_fusion_partner_segments {
    my ($cds_obj, $breakpoint_info) = @_;
    
    my ($chr, $breakpoint_coord, $orient) = split(/:/, $breakpoint_info);
    
    my ($left_segs_aref, $right_segs_aref) = &split_cds_at_breakpoint($cds_obj, $breakpoint_coord);

    if ($orient eq '+') {
        return(@$left_segs_aref);
    }
    else {
        return(reverse @$right_segs_aref);
    }
}

####
sub get_right_fusion_partner_segments {
    my ($cds_obj, $breakpoint_info) = @_;
    
    my ($chr, $breakpoint_coord, $orient) = split(/:/, $breakpoint_info);
    
    my ($left_segs_aref, $right_segs_aref) = &split_cds_at_breakpoint($cds_obj, $breakpoint_coord);
    
    if ($orient eq '+') {
        return(@$right_segs_aref);
    }
    else {
        return(reverse @$left_segs_aref);
    }
}



####
sub split_cds_at_breakpoint {
    my ($cds_obj, $breakpoint_coord) = @_;
    
    my @segments = sort {$a->{lend}<=>$b->{lend}} @{$cds_obj->{phased_segments}};
    
    my @segs_left;
    my @segs_right;

    
    foreach my $segment (@segments) {
        if ($segment->{rend} <= $breakpoint_coord) {
            push (@segs_left, $segment);
        }
        elsif ($segment->{lend} >= $breakpoint_coord) {
            push (@segs_right, $segment);
        }
        elsif(&overlaps_breakpoint($breakpoint_coord ,$segment->{lend}, $segment->{rend})) {

            ## split the segment at the breakpoint, keep breakpoint coordinate in each piece.
            my $orient = $segment->{orient};
            if ($orient eq '+') {
                my $new_left_segment = { chr => $segment->{chr}, 
                                         lend => $segment->{lend},
                                         rend => $breakpoint_coord,
                                         orient => $orient,
                                         rel_lend => $segment->{rel_lend},
                                         rel_rend => $segment->{rel_lend} + ($breakpoint_coord - $segment->{lend}),
                                         phase_beg => $segment->{phase_beg},
                                         phase_end => ($segment->{rel_lend} + $segment->{phase_beg} + ($breakpoint_coord - $segment->{rel_lend})) % 3,
                };
                
                my $new_right_segment = { chr => $segment->{chr},
                                          lend => $breakpoint_coord,
                                          rend => $segment->{rend},
                                          orient => $orient,
                                          rel_lend => $new_left_segment->{rel_rend},
                                          rel_rend => $segment->{rel_rend},
                                          phase_beg => $new_left_segment->{phase_end},
                                          phase_end => $segment->{phase_end},
                };

                push (@segs_left, $new_left_segment);
                push (@segs_right, $new_right_segment);
            }
            else {
                ## orient eq '-'
                
                my $new_right_segment = { chr => $segment->{chr},
                                          lend => $breakpoint_coord,
                                          rend => $segment->{rend},
                                          rel_lend => $segment->{rel_rend} + ($segment->{rend} - $breakpoint_coord),
                                          rel_rend => $segment->{rel_rend},
                                          phase_beg => $segment->{phase_beg},
                                          phase_end => ($segment->{rel_lend} + $segment->{phase_beg} + (($segment->{rend} - $breakpoint_coord)) % 3),
                };

                my $new_left_segment = { chr => $segment->{chr},
                                         lend => $segment->{lend},
                                         rend => $breakpoint_coord,
                                         rel_lend => $segment->{rel_lend},
                                         rel_rend => $new_right_segment->{rel_lend},
                                         phase_beg => $new_right_segment->{phase_end},
                                         phase_end => $segment->{phase_end},
                };
                
                push (@segs_left, $new_left_segment);
                push (@segs_right, $new_right_segment);
                
            }                
            
        }
        else {
            die "Error, shouldn't get here";
        }
    }
    
    return(\@segs_left, \@segs_right);
    
}

####
sub overlaps_breakpoint {
    my ($breakpoint_coord, $lend, $rend) = @_;

    if ($breakpoint_coord >= $lend && $breakpoint_coord <= $rend) {
        return(1);
    }
    else {
        return(0);
    }
}

####
sub segments_to_string {
    my (@segments) = @_;

    @segments = sort {$a->{lend}<=>$b->{lend}} @segments;

    my $chr = $segments[0]->{chr};
    my $orient = $segments[0]->{orient};

    my @coord_text;
    foreach my $segment (@segments) {
        push (@coord_text, join("-", $segment->{lend}, $segment->{rend}));
    }

    my $descr_text = join("|", $chr, $orient, @coord_text);

    return($descr_text);
    
}
