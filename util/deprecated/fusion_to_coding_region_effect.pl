#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Carp::Assert;
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);
use Data::Dumper;
use JSON::XS;
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use TiedHash;
use Nuc_translator;

my $DEBUG = 0;

my $usage = <<__EOUSAGE__;

##########################################################################
#
#  --fusions <string>                : fusion predictions
#
#  --genome_lib_dir <string>         : CTAT genome lib
#                                      
#  optional:
#
#  --show_all                        : by default, only the single 'best' fusion is shown
#                                      prioritizizing INFRAME and longer fusion transcripts.
#
##########################################################################


__EOUSAGE__

    ;


my $help_flag;
my $fusions_file;
my $genome_lib_dir;
my $show_all_flag = 0;

&GetOptions ( 'h' => \$help_flag,
              'fusions=s' => \$fusions_file,
              'genome_lib_dir=s' => \$genome_lib_dir,
   
              'show_all' => \$show_all_flag,
    );

unless ($fusions_file && $genome_lib_dir) {
    
    die $usage;
}

my $prot_info_db = "$genome_lib_dir/ref_annot.prot_info.dbm";

unless (-s $prot_info_db) {
    confess "Error, cannot locate $prot_info_db - be sure to use the latest version of the CTAT genome lib";
}


my %priority = ( 'INFRAME' => 1,
                 'FRAMESHIFT' => 2,
                 'NA' => 3,
                 '.' => 4);



main: {

    open (my $fh, $fusions_file) or die "Error, cannot open file $fusions_file";
    my $header = <$fh>;
    
    unless ($header =~ /^\#/) {
        die "Error, fusion_file: $fusions_file has unrecognizable header: $header";
    }
    chomp $header;
   
    my %header_to_index;
    {
        my @fields = split(/\t/, $header);
        for (my $i = 0; $i <= $#fields; $i++) {
            $header_to_index{$fields[$i]} = $i;
        }
        
        # ensure we have the fields we need:
        my @fields_need = ('#FusionName', 'LeftGene', 'LeftBreakpoint', 'RightGene', 'RightBreakpoint');
        my $missing_field = 0;
        foreach my $field (@fields_need) {
            unless (exists $header_to_index{$field}) {
                print STDERR "ERROR, missing column header: $field\n";
                $missing_field = 1;
            }
        }
        if ($missing_field) {
            die "Error, at least one column header was missing";
        }
    }
    

    print $header . "\t" . join("\t", 
                                "CDS_LEFT_ID", 
                                "CDS_LEFT_RANGE",
                                "CDS_RIGHT_ID",
                                "CDS_RIGHT_RANGE",
                                "PROT_FUSION_TYPE",
                                "FUSION_MODEL",
                                "FUSION_CDS",
                                "FUSION_TRANSL",
                                "PFAM_LEFT",
                                "PFAM_RIGHT",
        ) . "\n";
        
    
    my $prot_info_db_tied_hash = new TiedHash( { 'use' => $prot_info_db } );

    my $pfam_domain_dbm = "$genome_lib_dir/pfam_domains.dbm";
    my $pfam_domain_db_tied_hash = undef;
    if (-s $pfam_domain_dbm) {
        $pfam_domain_db_tied_hash = new TiedHash( { 'use' => $pfam_domain_dbm } );
    }
    
    
    while (<$fh>) {
        chomp;        
        my $line = $_;

        my @x = split(/\t/);

        my $fusion_name = $x[ $header_to_index{'#FusionName'} ];
        my $gene_left = $x[ $header_to_index{'LeftGene'}];
        my $break_left = $x[ $header_to_index{'LeftBreakpoint'} ];
                
        my $gene_right = $x[ $header_to_index{'RightGene'} ];
        my $break_right = $x[ $header_to_index{'RightBreakpoint'} ];
        

        # remove gene symbol, just want gene ID
        $gene_left =~ s/^.*\^//;
        $gene_right =~ s/^.*\^//;
        

        my @results = &examine_fusion_coding_effect($gene_left, $break_left, $gene_right, $break_right, $prot_info_db_tied_hash, $pfam_domain_db_tied_hash);
        
        if (@results) {
            ## just take the single 'best' one, arbitrarily chosen as the one with the longest fusion sequence.
            @results = sort { 
                $priority{$a->{prot_fusion_type}} <=> $priority{$b->{prot_fusion_type}}
                ||
                length($b->{prot_fusion_seq}) <=> length($a->{prot_fusion_seq}) } @results;
            
            #my $result = shift @results;
            
            foreach my $result (@results) {
                print join("\t", 
                           $line,
                           $result->{cds_left_id}, $result->{cds_left_range},
                           $result->{cds_right_id}, $result->{cds_right_range},
                           $result->{prot_fusion_type},
                           $result->{fusion_coding_descr},
                           $result->{cds_fusion_seq},
                           $result->{prot_fusion_seq},
                           
                           $result->{left_domains},
                           $result->{right_domains}
                           
                    ) . "\n";
                
                unless ($show_all_flag) {
                    # only showing first entry.
                    last;
                }
            }
        }
        else {
            print "$line" . ("\t." x 10) . "\n";
        }

    }
    close $fh;
    
    
    exit(0);
}


####
sub examine_fusion_coding_effect {
    my ($gene_left, $break_left, $gene_right, $break_right, $prot_info_db_tied_hash, $pfam_domain_db_tied_hash) = @_;
    
    my $gene_left_aref = &decode_my_json($gene_left, $prot_info_db_tied_hash, 1);
    my $gene_right_aref = &decode_my_json($gene_right, $prot_info_db_tied_hash, 1);
    
    unless ($gene_left_aref && $gene_right_aref) {
        return();
    }

    if (0) {
        ## debugging
        foreach my $cds_obj (@$gene_left_aref, @$gene_right_aref) {
            my @segments = @{$cds_obj->{phased_segments}};
            print &segments_to_string(@segments) . "\n";
        }
        
        return();
    }
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

                my $prot_fusion_type = "NA";
                if ($left_end_phase ne '.' && $right_beg_phase ne '.') {
                    $prot_fusion_type = ( ($left_end_phase + 1) % 3 == $right_beg_phase) ? "INFRAME" : "FRAMESHIFT";
                }
                
                my $left_segs_string = &segments_to_string(@$left_fuse_segments_aref);
                my $right_segs_string = &segments_to_string(@$right_fuse_segments_aref);

                my $left_domains_string = &get_pfam_domains($cds_left_id, $left_rel_rend, "left", $pfam_domain_db_tied_hash) || ".";
                my $right_domains_string = &get_pfam_domains($cds_right_id, $right_rel_lend, "right", $pfam_domain_db_tied_hash) || ".";
                
                
                push (@results, { cds_left_id => $cds_left_id,
                                  cds_right_id => $cds_right_id,
                                  cds_left_range => "1-$left_rel_rend",
                                  cds_right_range => "$right_rel_lend-" . length($cds_right_seq),
                                  prot_fusion_type => $prot_fusion_type,
                                  cds_fusion_seq => $fusion_seq,
                                  prot_fusion_seq => $pep,
                                  fusion_coding_descr => join("<==>", $left_segs_string, $right_segs_string),
                                  left_domains => $left_domains_string,
                                  right_domains => $right_domains_string,
                                  
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
    
    # ensure breakpoint overlaps a coding segment
    unless (&breakpoint_overlaps_cds_segment($cds_obj, $breakpoint_coord)) {
        return();
    }
    
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

    # ensure breakpoint overlaps a coding segment
    unless (&breakpoint_overlaps_cds_segment($cds_obj, $breakpoint_coord)) {
        return();
    }
    
    my ($left_segs_aref, $right_segs_aref) = &split_cds_at_breakpoint($cds_obj, $breakpoint_coord);
    
    if ($orient eq '+') {
        return(@$right_segs_aref);
    }
    else {
        return(reverse @$left_segs_aref);
    }
}



####
sub breakpoint_overlaps_cds_segment {
    my ($cds_obj, $breakpoint_coord) = @_;

    my @segments = sort {$a->{lend}<=>$b->{lend}} @{$cds_obj->{phased_segments}};

    foreach my $segment (@segments) {
        if ($segment->{lend} <= $breakpoint_coord && $breakpoint_coord <= $segment->{rend}) {
            return(1);
        }
    }

    return(0); # no overlap
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
                                         phase_end => ".", # set below
                };
                if ($segment->{phase_beg} ne '.') {
                    $new_left_segment->{phase_end} = ($segment->{phase_beg} + ($breakpoint_coord - $segment->{lend})) % 3;
                }
                
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
                
                
                print "Original segment: " . &segments_to_string($segment) . 
                    "\nNew frags: " . &segments_to_string($new_left_segment) . "\t" . &segments_to_string($new_right_segment) . "\n\n" if $DEBUG;
                
                
                if ($segment->{phase_beg} ne '.') {
                    my $supposed_end_phase = ($new_right_segment->{phase_beg} + ($new_right_segment->{rel_rend} - $new_right_segment->{rel_lend})) % 3;
                    print "supposed end phase: $supposed_end_phase\n" if $DEBUG;
                    assert($supposed_end_phase == $new_right_segment->{phase_end});
                }
                
             
            }
            else {
                ## orient eq '-'
                
                my $new_right_segment = { chr => $segment->{chr},
                                          lend => $breakpoint_coord,
                                          rend => $segment->{rend},
                                          orient => $orient,
                                          rel_lend => $segment->{rel_rend} + ($segment->{rend} - $breakpoint_coord),
                                          rel_rend => $segment->{rel_rend},
                                          phase_beg => $segment->{phase_beg},
                                          phase_end => ".", # set below
                };
                
                if ($segment->{phase_beg} ne '.') {
                    $new_right_segment->{phase_end} = ($segment->{phase_beg} + (($segment->{rend} - $breakpoint_coord))) % 3;
                }

                my $new_left_segment = { chr => $segment->{chr},
                                         lend => $segment->{lend},
                                         rend => $breakpoint_coord,
                                         orient => $orient,
                                         rel_lend => $segment->{rel_lend},
                                         rel_rend => $new_right_segment->{rel_lend},
                                         phase_beg => $new_right_segment->{phase_end},
                                         phase_end => $segment->{phase_end},
                };
                
                push (@segs_left, $new_left_segment);
                push (@segs_right, $new_right_segment);
                
                print "Original segment: " . &segments_to_string($segment) . 
                    "\nNew frags: " . &segments_to_string($new_left_segment) . "\t" . &segments_to_string($new_right_segment) . "\n\n" if $DEBUG;
                
                if ($segment->{phase_beg} ne ".") {
                    assert($new_right_segment->{phase_end} == ($new_right_segment->{phase_beg} + $new_right_segment->{rend} - $new_right_segment->{lend}) % 3);
                    
                    assert($new_left_segment->{phase_end} == ($new_left_segment->{phase_beg} + $new_left_segment->{rend} - $new_left_segment->{lend}) % 3);
                }
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
    my $orient = $segments[0]->{orient} or die "Error: " . Dumper(\@segments);

    my @coord_text;
    foreach my $segment (@segments) {
        my ($phase_left, $phase_right) = ($orient eq '+') 
            ? ($segment->{phase_beg}, $segment->{phase_end})
            : ($segment->{phase_end}, $segment->{phase_beg});
        
        push (@coord_text, "[$phase_left]" . join("-", $segment->{lend}, $segment->{rend}) . "[$phase_right]");
    }
    
    my $descr_text = join("|", $chr, $orient, @coord_text);

    return($descr_text);
    
}


####
sub get_pfam_domains {
    my ($cds_id, $cds_coord, $left_or_right_side, $pfam_domain_db_tied_hash) = @_;

    unless ($pfam_domain_db_tied_hash) {
        return;
    }

    my $pfam_hits_aref = &decode_my_json($cds_id, $pfam_domain_db_tied_hash);

    unless ($pfam_hits_aref) { return; }
    
    my @pfam_hits = @$pfam_hits_aref;
    
    my @pfam_domains_selected;

    foreach my $pfam_hit (@pfam_hits) {
        my ($start, $end) = ($pfam_hit->{query_start}, $pfam_hit->{query_end});

        if  (($left_or_right_side eq 'left' && $end <= $cds_coord)
             ||
             ($left_or_right_side eq 'right' && $start >= $cds_coord) ) {

            ## domain entirely on the side of the protein included in the fusion.
            
            push (@pfam_domains_selected, $pfam_hit);
        }
        elsif ($start < $cds_coord && $cds_coord < $end) {
            ## overlaps

            ## fragment it and return the fragment.
            my $pfam_hit_copy = &clone($pfam_hit);
            if ($left_or_right_side eq 'left') {
                $pfam_hit_copy->{query_end} = $cds_coord;
                $pfam_hit_copy->{query_end_partial} = 1;
            }
            else {
                # right side
                $pfam_hit_copy->{query_start} = $cds_coord;
                $pfam_hit_copy->{query_start_partial} = 1;
            }
            $pfam_hit_copy->{hmmer_domain} .= "-PARTIAL";
            push (@pfam_domains_selected, $pfam_hit_copy);
            
        }
    }

    ## generate a summary string.
    @pfam_domains_selected = sort {$a->{query_start}<=>$b->{query_start}} @pfam_domains_selected;

    my @pfam_descrs;
    foreach my $pfam_domain (@pfam_domains_selected) {
        #'query_start' => '369',
        #'cds_id' => 'DISP1|ENST00000284476.6',
        #'domain_evalue' => '7.3e-21',
        #'query_end' => '733',
        #'hmmer_domain' => 'Patched'

        if ($pfam_domain->{query_start_partial}) {
            $pfam_domain->{query_start} = "~" . $pfam_domain->{query_start};
        }
        if ($pfam_domain->{query_end_partial}) {
            $pfam_domain->{query_end} = $pfam_domain->{query_end} . "~";
        }
        
        my $descr = join("|", $pfam_domain->{hmmer_domain},
                         $pfam_domain->{query_start} . "-" . $pfam_domain->{query_end},
                         $pfam_domain->{domain_evalue});

        push (@pfam_descrs, $descr);
    }

    my $ret_descr = join("^", @pfam_descrs);

    return($ret_descr);
            
}


####
sub clone {
    my ($hashref) = @_;

    my $clone_ref = {};
    foreach my $key (keys %$hashref) {

        $clone_ref->{$key} = $hashref->{$key};
    }

    return($clone_ref);
}


####
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

