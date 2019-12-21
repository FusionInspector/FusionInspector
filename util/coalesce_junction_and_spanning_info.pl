#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use Carp;
use lib ("$FindBin::Bin/../PerlLib");
use DelimParser;
use Data::Dumper;

my $PSEUDOCOUNT = 0;

my $usage = "\n\tusage: $0 junction_info_A.txt,[junction_info_B.txt,...] spanning_info_A.txt,[spanning_info_B.txt,...] [PSEUDOCOUNT=$PSEUDOCOUNT]\n\n";

my $junction_info_file_list = $ARGV[0] or die $usage;
my $spanning_info_file_list = $ARGV[1] or die $usage;
if ($ARGV[2]) {
    $PSEUDOCOUNT = $ARGV[2];
}


main: {

    my %fusion_info;
    my %fusion_large_breakpoint_anchored;

    foreach my $junction_info_file (split(/,/, $junction_info_file_list)) {

        &parse_junction_info_file(\%fusion_info, "junction", $junction_info_file, \%fusion_large_breakpoint_anchored);
    }
    
    foreach my $spanning_info_file (split(/,/, $spanning_info_file_list)) {

        &parse_spanning_info_file(\%fusion_info, "spanning", $spanning_info_file);

    }


    my @fields = ("#FusionName", "JunctionReadCount", "SpanningFragCount", 
                  "LeftGene", "LeftLocalBreakpoint", "LeftBreakpoint",
                  "RightGene", "RightLocalBreakpoint", "RightBreakpoint",
                  "SpliceType", "LargeAnchorSupport", 
                  "JunctionReads", "SpanningFrags",
                  "NumCounterFusionLeft", "CounterFusionLeftReads",
                  "NumCounterFusionRight", "CounterFusionRightReads",
                  "FAR_left", "FAR_right");
    
    my $tab_writer = new DelimParser::Writer(*STDOUT, "\t", \@fields);

    foreach my $fusion (keys %fusion_info) {

        
        
        my @junction_reads = keys %{$fusion_info{$fusion}->{'junction'}};

        my $has_large_anchor_junction_support = $fusion_large_breakpoint_anchored{$fusion} || "NO";
        
        my @spanning_reads;
        if (exists $fusion_info{$fusion}->{'spanning'}) {
            @spanning_reads = keys %{$fusion_info{$fusion}->{'spanning'}};
        
            @spanning_reads = &remove_residual_junction_from_spanning_reads(\@junction_reads, \@spanning_reads);
            
        }
        
        my $num_junction_reads = &count_frags(@junction_reads); # need to be careful when read pairs overlap and they're both split reads.
        my $num_spanning_reads = scalar(@spanning_reads);
        
        my @left_contrary_reads;
        if (exists $fusion_info{$fusion}->{'spanning-left_contrary_support'}) {
            @left_contrary_reads = keys %{$fusion_info{$fusion}->{'spanning-left_contrary_support'}};
        }
        my $num_left_contrary_reads = scalar(@left_contrary_reads);
        if ($num_left_contrary_reads == 0) {
            push (@left_contrary_reads, ".");
        }
        

        my @right_contrary_reads;
        if (exists $fusion_info{$fusion}->{'spanning-right_contrary_support'}) {
            @right_contrary_reads = keys %{$fusion_info{$fusion}->{'spanning-right_contrary_support'}};
        }
        my $num_right_contrary_reads = scalar(@right_contrary_reads);
        if ($num_right_contrary_reads == 0) {
            push (@right_contrary_reads, ".");
        }
        
        my $FAR_left = ($num_junction_reads + $num_spanning_reads + $PSEUDOCOUNT) 
            / 
            ($num_left_contrary_reads + $PSEUDOCOUNT);
        
        $FAR_left = sprintf("%.2f", $FAR_left);
        
        my $FAR_right = ($num_junction_reads + $num_spanning_reads + $PSEUDOCOUNT) 
            / 
            ($num_right_contrary_reads + $PSEUDOCOUNT);
        
        $FAR_right = sprintf("%.2f", $FAR_right);
        
        my ($geneA, $local_brkpt_A, $chr_brkpt_A, $geneB, $local_brkpt_B, $chr_brkpt_B, $spliceType) = split(/\t/, $fusion);
            
        my ($geneA_symbol, @restA) = split(/\^/, $geneA);
        my ($geneB_symbol, @restB) = split(/\^/, $geneB);
        
        $tab_writer->write_row( { 

            '#FusionName' => "$geneA_symbol--$geneB_symbol",
            LeftGene => $geneA,
            LeftLocalBreakpoint => $local_brkpt_A,
            LeftBreakpoint => $chr_brkpt_A,
            RightGene => $geneB,
            RightLocalBreakpoint => $local_brkpt_B,
            RightBreakpoint => $chr_brkpt_B,
            SpliceType => $spliceType,
            JunctionReadCount => $num_junction_reads,
            SpanningFragCount => $num_spanning_reads,
            LargeAnchorSupport => $has_large_anchor_junction_support,
            JunctionReads => join(",", @junction_reads),
            SpanningFrags => join(",", @spanning_reads),
            NumCounterFusionLeft => $num_left_contrary_reads,
            CounterFusionLeftReads => join(",", @left_contrary_reads),
            NumCounterFusionRight => $num_right_contrary_reads,
            CounterFusionRightReads => join(",", @right_contrary_reads),
            FAR_left => $FAR_left,
            FAR_right => $FAR_right,
                                } );
        
    }
    
    

    exit(0);
}


####
sub parse_junction_info_file {
    my ($fusion_info_href, $fusion_read_type, $file, $fusion_large_breakpoint_anchored_href) = @_;

    open (my $fh, $file) or die "Error, cannot open file $file";

    my $tab_reader = new DelimParser::Reader($fh, "\t");
    
    while (my $row = $tab_reader->get_row()) {
        my ($geneA, $coordA, $orig_coordA, 
            $geneB, $coordB, $orig_coordB, 
            $splice_info, $count, $has_large_anchors, $read_list) = ($row->{LeftGene},
                                                                     $row->{LeftLocalBreakpoint},
                                                                     $row->{LeftBreakpoint},
                                                                     $row->{RightGene},
                                                                     $row->{RightLocalBreakpoint},
                                                                     $row->{RightBreakpoint},
                                                                     $row->{SpliceType},
                                                                     $row->{JunctionReadCount},
                                                                     $row->{LargeAnchorSupport},
                                                                     $row->{JunctionReads});
        
        
        my $fusion_token = join("\t", $geneA, $coordA, $orig_coordA, $geneB, $coordB, $orig_coordB, $splice_info);
        
        foreach my $read (split(/,/, $read_list)) {

            $fusion_info_href->{$fusion_token}->{$fusion_read_type}->{$read}++;
        }
        
        $fusion_large_breakpoint_anchored_href->{$fusion_token} = $has_large_anchors;
    }
    close $fh;
    
    return;
}



####
sub parse_spanning_info_file {
    my ($fusion_info_href, $fusion_read_type, $file) = @_;

    open (my $fh, $file) or die "Error, cannot open file $file";
    
    my $tab_reader = new DelimParser::Reader($fh, "\t");

    while (my $row = $tab_reader->get_row()) {
        
        my ($geneA, $coordA, $orig_coordA, 
            $geneB, $coordB, $orig_coordB, 
            $splice_info, $count, $read_list,
            $left_contrary_support_count, $left_contrary_support_reads,
            $right_contrary_support_count, $right_contrary_support_reads) = ($row->{LeftGene},
                                                                             $row->{LeftLocalBreakpoint},
                                                                             $row->{LeftBreakpoint},
                                                                             $row->{RightGene},
                                                                             $row->{RightLocalBreakpoint},
                                                                             $row->{RightBreakpoint},
                                                                             $row->{SpliceType},
                                                                             $row->{SpanningFragCount},
                                                                             $row->{SpanningFrags},
                                                                             $row->{NumCounterFusionLeft},
                                                                             $row->{CounterFusionLeftReads},
                                                                             $row->{NumCounterFusionRight},
                                                                             $row->{CounterFusionRightReads});
        
        my $fusion_token = join("\t", $geneA, $coordA, $orig_coordA, $geneB, $coordB, $orig_coordB, $splice_info);
        
        foreach my $read (split(/,/, $read_list)) {

            $fusion_info_href->{$fusion_token}->{$fusion_read_type}->{$read}++;
        }


        if ($left_contrary_support_count > 0) {
            foreach my $read (split(/,/, $left_contrary_support_reads)) {
                
                $fusion_info_href->{$fusion_token}->{"$fusion_read_type-left_contrary_support"}->{$read} = 1;
            }
        }
            
        if ($right_contrary_support_count > 0) {
            foreach my $read (split(/,/, $right_contrary_support_reads)) {
                
                $fusion_info_href->{$fusion_token}->{"$fusion_read_type-right_contrary_support"}->{$read} = 1;
                
            }
        }
    }
    
    close $fh;

    
    return;
}


####
sub remove_residual_junction_from_spanning_reads {
    my ($junction_reads_aref, $spanning_reads_aref) = @_;

    # junction takes precendence over spanning.
    
    my %junction;
    foreach my $junction_read (@$junction_reads_aref) {
        
        my $read_to_discard = $junction_read;
        $read_to_discard =~ s|/[12]$||;
        
        $junction{$read_to_discard}++;
    }

    my @spanning_retain;
    
    foreach my $spanning_read (@$spanning_reads_aref) {
        if (exists $junction{$spanning_read}) {
            print STDERR "-discarding spanning read: $spanning_read since found among merged junction reads.\n";
        }
        else {
            push (@spanning_retain, $spanning_read);
        }
    }

    return(@spanning_retain);
}

####
sub count_frags {
    my @read_names = @_;
    
    my %frags;

    foreach my $read_name (@read_names) {
        $read_name =~ s/\/[12]$//;
        $frags{$read_name}++;
    }

    my $count = scalar(keys %frags);

    return($count);
}


