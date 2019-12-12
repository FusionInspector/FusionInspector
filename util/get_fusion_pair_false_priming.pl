#!/usr/bin/env perl

use strict;
use warnings;
use Carp;

use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Process_cmd;
use Fasta_reader;
use DelimParser;
use Data::Dumper;
use Set::IntervalTree;
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);


my $min_reputer_len = 11;

my $TMPDIR = $ENV{TMPDIR} || "/tmp";

my $BREAKPOINT_RANGE = 25;

my $usage = <<__EOUSAGE__;

###############################################################
#
# Required:
#
#  --gtf_file <string>              finspector.gtf
#  --fasta_file <string>            finspector.fa
#  --junction_info_file <string>    file.fusion_junction_info 
#
# Optional:
#
#  --output <string>           output filename (default: writes to stdout)
#
#  --min_reputer_len <int>     default: $min_reputer_len
#  --breakpoint_range <int>    distance allowed around breakpoints (default: $BREAKPOINT_RANGE)
#
#  --DEBUG|d
#  --TMPDIR <string>           default: $TMPDIR
#
##############################################################



__EOUSAGE__

    ;

my $help_flag;
my $DEBUG = 0;
my $gtf_file;
my $fasta_file;
my $output;
my $junction_info_file;

&GetOptions('help|h' => \$help_flag,
            
            'gtf_file=s' => \$gtf_file,
            'fasta_file=s' => \$fasta_file,
            'junction_info_file=s' => \$junction_info_file,
            
            'output=s' => \$output, 
            
            'min_reputer_len=i' => \$min_reputer_len,
            'breakpoint_range=i' => \$BREAKPOINT_RANGE,
            
            'DEBUG|d' => \$DEBUG,
            'TMPDIR=s' => \$TMPDIR,
    );

if ($help_flag) {
    die $usage;
}

unless ($gtf_file && $fasta_file && $junction_info_file) {
    die $usage;
}

my $ofh;
if ($output) {
    open($ofh, ">$output") or die "Error, cannot write to file: $output";
}
else {
    $ofh = *STDOUT;
}



main: {


    ##############################################################
    ## Get the reference gene coordinates on each fusion-scaffold

    my %scaffold_itree_to_exon_structs;
    
    my %scaffold_to_gene_structs = &parse_gtf_file($gtf_file, \%scaffold_itree_to_exon_structs);
    
    
    my %scaffold_to_gene_inner_bounds;
    
    {
        foreach my $scaffold (keys %scaffold_to_gene_structs) {

            my @gene_structs = @{$scaffold_to_gene_structs{$scaffold}};
            
            @gene_structs = sort {$a->{lend} <=> $b->{lend}} @gene_structs;

            if (scalar @gene_structs > 2) {
                @gene_structs = &partition_gene_structs(@gene_structs);
            }
                        
            
            my $left_gene = $gene_structs[0];
            my $right_gene = $gene_structs[1];
            
            my ($geneA, $geneB) = ($left_gene->{gene_id}, $right_gene->{gene_id});
            
            my $gene_bound_left = $left_gene->{rend};
            my $gene_bound_right = $right_gene->{lend};
            
            
            if ($gene_bound_left > $gene_bound_right) {
                die "Error, bounds out of order: $gene_bound_left ! <  $gene_bound_right";
            }

            $scaffold_to_gene_inner_bounds{$scaffold} = [$gene_bound_left, $gene_bound_right];
            
        }
    }

    my %scaffold_to_junction_breakpoints;
    {
        open(my $fh, $junction_info_file) or die "Error, cannot open file: $junction_info_file";
        my $delim_reader = new DelimParser::Reader($fh, "\t");
        while(my $row = $delim_reader->get_row()) {

            my $left_gene_info = $delim_reader->get_row_val($row, "LeftGene");
            my @left_gene_pts = split('\^', $left_gene_info);
            my $left_gene = $left_gene_pts[0];

            my $right_gene_info = $delim_reader->get_row_val($row, "RightGene");
            my @right_gene_pts = split('\^', $right_gene_info);
            my $right_gene = $right_gene_pts[0];

            my $fusion_scaffold = "$left_gene--$right_gene";
            
            my $left_breakpoint = $delim_reader->get_row_val($row, "LeftLocalBreakpoint");
            my $right_breakpoint = $delim_reader->get_row_val($row, "RightLocalBreakpoint");

            push (@{$scaffold_to_junction_breakpoints{$fusion_scaffold}}, [$left_breakpoint, $right_breakpoint]);
        }
    }
    
    
    ## run reputer
    my $fasta_reader = new Fasta_reader($fasta_file);
    while (my $seq_obj = $fasta_reader->next()) {
        my $contig_acc = $seq_obj->get_accession();
        
        my $gene_bounds_aref = $scaffold_to_gene_inner_bounds{$contig_acc} or die "Error, no bounds defined for $contig_acc";
        my ($bound_left, $bound_right) = @$gene_bounds_aref;

        my $itree = $scaffold_itree_to_exon_structs{$contig_acc} or die "Error, no interval tree for $contig_acc";

        my $fusion_breakpoints_aref = $scaffold_to_junction_breakpoints{$contig_acc};
        unless ($fusion_breakpoints_aref) {
            # no fusion calls here... skip
            next;
        }
        
        
        my $fasta_entry = $seq_obj->get_FASTA_format();
        my $seq_filename = "$$.seq.fa";
        my $tmp_fasta_file = "$TMPDIR/$seq_filename";
        {
            open(my $ofh, ">$tmp_fasta_file") or confess("Error, cannot write to file: $tmp_fasta_file");
            print $ofh $fasta_entry;
            close $ofh;
        }
        
        my $cmd = "mkvtree -db $tmp_fasta_file -dna -pl -allout -v >/dev/null 2>&1";
        &process_cmd($cmd);

        $cmd = "mv $$.seq.fa.* $TMPDIR/";
        &process_cmd($cmd);
        
        my $vmatch_outfile = "$TMPDIR/$seq_filename.vmatch.out";
        $cmd = "vmatch -v -l $min_reputer_len -e 1 -showdesc 10  $tmp_fasta_file > $vmatch_outfile 2>/dev/null";
        &process_cmd($cmd);
        
        open(my $fh, $vmatch_outfile) or die "Error, cannot read file: $vmatch_outfile";
        my $match_counter = 0;
        while(<$fh>) {
            if (/^\#/) {
                next;
            }
            chomp;
            s/^\s+//;
            my @x = split(/\s+/);
            my $geneA_match_lend = $x[2];
            my $geneA_match_length = $x[0];
            my $geneA_match_rend = $geneA_match_lend + $geneA_match_length - 1;

            my $geneB_match_lend = $x[6];
            my $geneB_match_length = $x[4];
            my $geneB_match_rend = $geneB_match_lend + $geneB_match_length - 1;

            # only want the matches between the fusion gene pairs
            unless ($geneA_match_lend < $bound_left && $geneB_match_rend > $bound_right) { next; }

            # make sure it overlaps exons for each gene
            unless (&overlaps_an_exon($geneA_match_lend, $geneA_match_rend, $itree)
                    &&
                    &overlaps_an_exon($geneB_match_lend, $geneB_match_rend, $itree) ) {

                next;
            }

            unless (&near_fusion_breakpoint([$geneA_match_lend, $geneA_match_rend],
                                            [$geneB_match_lend, $geneB_match_rend],
                                            $fusion_breakpoints_aref,
                                            $BREAKPOINT_RANGE)) {
                next;
            }
            
            # output bed entry
                        
            print $ofh join("\t", 
                            $contig_acc,
                            $geneA_match_lend - 1, # zero based coord
                            $geneB_match_rend, # rend is exclusive
                            "match_" . ++$match_counter . "_len_" . $geneA_match_length,
                            0,
                            "+",
                            $geneA_match_lend - 1,
                            $geneB_match_rend,
                            ".",
                            2,
                            join(",", $geneA_match_length, $geneB_match_length),
                            join(",", 0, $geneB_match_lend - $geneA_match_lend),
                ) . "\n";
            
        }

        close $fh;
        unlink (<$TMPDIR/$seq_filename*>);
        
        
        
    }
    
    close $ofh;

    exit(0);
}




sub parse_gtf_file {
    my ($gtf_file, $scaffold_itree_to_exon_structs_href) = @_;
    
    my %scaff_to_gene_to_coords;

    open (my $fh, $gtf_file) or die "Error, cannot open file $gtf_file";
    while (<$fh>) {
        chomp;
        my @x = split(/\t/);
        my $scaffold_id = $x[0];
        my $type = $x[2];
        
        unless ($type eq 'exon') { next; }
        
        my $interval_tree = $scaffold_itree_to_exon_structs_href->{$scaffold_id};
        unless (ref $interval_tree) {
            $interval_tree = $scaffold_itree_to_exon_structs_href->{$scaffold_id} = Set::IntervalTree->new;
        }
        
        my $info = $x[8];
        $info =~ /FI_gene_label \"([^\"]+)\"/ or die "Error, cannot parse FI_gene_label from $info";
        my $gene_id = $1;

        my ($lend, $rend) = ($x[3], $x[4]);
        push (@{$scaff_to_gene_to_coords{$scaffold_id}->{$gene_id}}, $lend, $rend);
     

        if ($lend != $rend) {
            $interval_tree->insert("$scaffold_id:$lend-$rend", $lend, $rend);
        }
        
    }
    close $fh;

    
    my %scaffold_to_gene_structs;

    foreach my $scaffold (keys %scaff_to_gene_to_coords) {
        my @genes = keys %{$scaff_to_gene_to_coords{$scaffold}};
    
        my @gene_structs;
    
        foreach my $gene (@genes) {
            my @coords = sort {$a<=>$b} @{$scaff_to_gene_to_coords{$scaffold}->{$gene}};
            my $lend = shift @coords;
            my $rend = pop @coords;
            push (@{$scaffold_to_gene_structs{$scaffold}}, { gene_id => $gene,
                                                             lend => $lend,
                                                             rend => $rend,
                  });
        }
        
    }
        
    return(%scaffold_to_gene_structs);
}

sub partition_gene_structs {
    my (@gene_structs) = @_;

    # should already be sorted by lend.

    my $left_gene_struct = shift @gene_structs;
    my @left_gene_structs = ($left_gene_struct);

    my $left_gene_id = $left_gene_struct->{gene_id};
    my $left_gene_symbol = $left_gene_id;
    $left_gene_symbol =~ s/\^.*$//;

    my @right_gene_structs;

    foreach my $gene_struct (@gene_structs) {
        my $gene_id = $gene_struct->{gene_id};
        my ($gene_sym, $rest) = split(/\^/, $gene_id);

        if ($gene_sym eq $left_gene_symbol) {
            push (@left_gene_structs, $gene_struct);
        }
        else {
            push (@right_gene_structs, $gene_struct);
        }
    }
 unless (@left_gene_structs && @right_gene_structs) {
        die "Error, couldn't partition gene structs across contig breakpoint: " . Dumper(\@gene_structs);
    }

    @left_gene_structs = sort {$a->{rend}<=>$b->{rend}} @left_gene_structs;
    @right_gene_structs = sort {$a->{lend}<=>$b->{lend}} @right_gene_structs;

    # want those that are most adjacent to the central fusion contig join
    my @gene_struct_tuple = ($left_gene_structs[$#left_gene_structs], $right_gene_structs[0]);

    return(@gene_struct_tuple);
}

####
sub overlaps_an_exon {
    my ($lend, $rend, $itree) = @_;

    my $hits_aref = $itree->fetch($lend, $rend);

    if ($hits_aref && @$hits_aref) {
        return(1); # yes overlaps
    }
    else {
        return(0); # negatory
    }
}

####
sub near_fusion_breakpoint {
    my ($repeat_match_A_aref, $repeat_match_B_aref, $fusion_breakpoints_aref, $breakpoint_range) = @_;


    my ($A_match_lend, $A_match_rend) = @$repeat_match_A_aref;
    my ($B_match_lend, $B_match_rend) = @$repeat_match_B_aref;

    foreach my $fusion_breakpoint_aref (@$fusion_breakpoints_aref) {
        my ($fusion_lend, $fusion_rend) = @$fusion_breakpoint_aref;

        if ($A_match_lend < $fusion_lend + $breakpoint_range
            &&
            $A_match_rend > $fusion_lend - $breakpoint_range

            &&

            $B_match_lend < $fusion_rend + $breakpoint_range
            &&
            $B_match_rend > $fusion_rend - $breakpoint_range) {

            return(1); # yes, near breakpoint
        }
    }


    return(0); # no, not near a breakpoint
}


