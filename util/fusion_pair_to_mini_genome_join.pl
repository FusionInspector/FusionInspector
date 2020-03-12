#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Long qw(:config no_ignore_case bundling pass_through);

use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Nuc_translator;
use Overlap_piler;
use Data::Dumper;

my $max_intron_length = 1000;
my $genome_flank_size = 1000;

my $usage = <<__EOUSAGE__;

###############################################################################################
#
#  Required:
#
#  --fusions <string>               file containing list of fusion pairs (format:  geneA--geneB)
#
#  --gtf <string>                   genome annotation in gtf format 
#
#  --genome_fa <string>             genome sequence in fasta format
#
# Optional:
#
#  --shrink_introns
#
#  --max_intron_length <int>        default: $max_intron_length  (only when --shrink_introns used)
#
#  --genome_flank <int>             amt. of genomic sequence to extract flanking each gene (default: $genome_flank_size)
#
#  --out_prefix <string>            output prefix for output files (gtf and fasta) default: geneMergeContig.\${process_id}
#
###############################################################################################


__EOUSAGE__

    ;

my $help_flag;

my $fusions_file;
my $gtf_file;
my $genome_fasta_file;
my $out_prefix = "geneMergeContig.$$";
my $shrink_introns_flag = 0;

&GetOptions ( 'h' => \$help_flag,
              
              'fusions=s' => \$fusions_file,
              
              'gtf=s' => \$gtf_file,
              'genome_fa=s' => \$genome_fasta_file,

              'shrink_introns' => \$shrink_introns_flag,
              'max_intron_length=i' => \$max_intron_length,
              'genome_flank=i' => \$genome_flank_size,
              
              'out_prefix=s' => \$out_prefix,
              

    );


if ($help_flag) {
    die $usage;
}

unless ($fusions_file && $gtf_file && $genome_fasta_file) {
    die $usage;
}

main: {

    my @chim_pairs;


    
  parse_fusion_candidates: {
      
      my @fusion_files = split(/,/, $fusions_file);
      foreach my $file (@fusion_files) {
          
          open (my $fh, $file) or die "Error, cannot open file $file";
          while (<$fh>) {
              if (/^\#/) { next; }
              unless (/\w/) { next; }
              
              my ($chim_pair, @rest) = split(/\s+/);
              if ($chim_pair =~ /^(\S+)--(\S+)$/) {
                  
                  my ($geneA, $geneB) = ($1, $2);
                  
                  if ($geneA eq $geneB) { next; } ## no selfies
                  
                  push (@chim_pairs, $chim_pair);
              }
              else {
                  die "Error, cannot parse $chim_pair as a fusion-gene candidate.";
              }
          }
          close $fh;
      }
    }
    
    my %genes_want;
    {
        foreach my $chim_pair (@chim_pairs) {
            foreach my $gene (split(/--/, $chim_pair)) {
                $genes_want{$gene} = 1;
                my @parts = split(/-/, $gene);
                if (scalar @parts > 1) {
                    foreach my $part (@parts) {
                        $genes_want{$part} = 1;
                    }
                }
            }
        }
    }
    
    
    my %gene_to_gtf = &extract_gene_gtfs($gtf_file, \%genes_want);
    
    
    ## split readthru transcripts into their separate parts
    my @tmp_chim_pairs;
    foreach my $chim_pair (@chim_pairs) {
        my ($left_gene, $right_gene) = split(/--/, $chim_pair);
        
        my @pairs;

        my @left_genes = split(/-/, $left_gene);
        my @right_genes = split(/-/, $right_gene);
        
        my $all_ok = 1;

        foreach my $tmp_left_gene (@left_genes) {
            
            foreach my $tmp_right_gene (@right_genes) {
                
                if (exists $gene_to_gtf{$tmp_left_gene}
                    &&
                    exists $gene_to_gtf{$tmp_right_gene}

                    &&
                    $tmp_left_gene ne $tmp_right_gene
                    ) {
                    
                    push (@pairs, "$tmp_left_gene--$tmp_right_gene");
                }
                else {
                    $all_ok = 0;
                }
            }
        }
        
        if ($all_ok && @tmp_chim_pairs) {
            push (@tmp_chim_pairs, @pairs);
        }
        else {
            # keep original version
            push (@tmp_chim_pairs, $chim_pair);
        }
    }
    
    @chim_pairs = @tmp_chim_pairs;
    

    open (my $out_genome_ofh, ">$out_prefix.fa.tmp") or die "Error, cannot write to $out_prefix.fa.tmp";
    open (my $out_gtf_ofh, ">$out_prefix.gtf.tmp") or die "Error, cannot write to $out_prefix.gtf.tmp";
    
    my %seen;
    
    
    foreach my $chim_pair (sort @chim_pairs) {
        
        if ($seen{$chim_pair}) {
            next; 
        }
        $seen{$chim_pair}++;
        
        my ($left_gene, $right_gene) = split(/--/, $chim_pair);
                        
        my $left_gene_gtf = $gene_to_gtf{$left_gene};
        my $right_gene_gtf = $gene_to_gtf{$right_gene};
        
        unless ($left_gene_gtf) {
            print STDERR "WARNING, no gtf annotations found for [$left_gene]\n";
            next;
        }
        unless ($right_gene_gtf) {
            print STDERR "WARNING, no gtf annotations found for [$right_gene]\n";
            next;
        }
        
        eval {

            my ($left_gene_supercontig_gtf, $left_gene_sequence_region) = &get_gene_contig_gtf($left_gene_gtf, $genome_fasta_file);
            
            my ($right_gene_supercontig_gtf, $right_gene_sequence_region) = &get_gene_contig_gtf($right_gene_gtf, $genome_fasta_file);
            
            if ($shrink_introns_flag) {
                ($left_gene_supercontig_gtf, $left_gene_sequence_region) = &shrink_introns($left_gene_supercontig_gtf, $left_gene_sequence_region, $max_intron_length);

                ($right_gene_supercontig_gtf, $right_gene_sequence_region) = &shrink_introns($right_gene_supercontig_gtf, $right_gene_sequence_region, $max_intron_length);
            }

            my $supercontig = $left_gene_sequence_region . ("N" x 1000) . $right_gene_sequence_region;
            
            $right_gene_supercontig_gtf = &adjust_gtf_coordinates($right_gene_supercontig_gtf, length($left_gene_sequence_region) + 1000);
            
            $supercontig =~ s/(\S{60})/$1\n/g; # make fasta 
            chomp $supercontig;
            
            print $out_genome_ofh ">$chim_pair\n$supercontig\n";
            
            
            my $out_gtf = $left_gene_supercontig_gtf . $right_gene_supercontig_gtf;
            $out_gtf = &set_gtf_scaffold_name($chim_pair, $out_gtf);
                        
            print $out_gtf_ofh $out_gtf;
            
        };

        if ($@) {
            print STDERR "$@\n";
        }
    }
    
    
    print STDERR "Done.\n";
    
    close $out_genome_ofh;
    close $out_gtf_ofh;

    if (! -s "$out_prefix.fa.tmp") {
        die "Error, no fusion contigs written";
    }
    else {
        rename("$out_prefix.fa.tmp", "$out_prefix.fa");
        rename("$out_prefix.gtf.tmp", "$out_prefix.gtf");
    }
    
    exit(0);
    

}

####
sub shrink_introns {
    my ($gene_gtf, $gene_seq_region, $max_intron_length) = @_;

    my @gtf_structs;
    my @gtf_lines = split(/\n/, $gene_gtf);
    foreach my $gtf_line (@gtf_lines) {
        my @x = split(/\t/, $gtf_line);
        push (@gtf_structs, \@x);
    }
    
    @gtf_structs = sort {$a->[3]<=>$b->[3]} @gtf_structs;

    
    ## get exon piles
    my $overlap_piler = new Overlap_piler();
    foreach my $gtf_row_aref (@gtf_structs) {
        
        my $exon_lend = $gtf_row_aref->[3];
        my $exon_rend = $gtf_row_aref->[4];
        $overlap_piler->add_coordSet($gtf_row_aref, $exon_lend, $exon_rend);
    }

    my @piles = $overlap_piler->build_clusters();
    
    my @pile_structs;
    foreach my $pile (@piles) {
        
        my @all_coords;
        foreach my $gtf_row_aref (@$pile) {
            my $lend = $gtf_row_aref->[3];
            my $rend = $gtf_row_aref->[4];
            push (@all_coords, $lend, $rend);
        }
        @all_coords = sort {$a<=>$b} @all_coords;

        my $pile_lend = shift @all_coords;
        my $pile_rend = pop @all_coords;
        
        my $pile_struct = { pile => $pile,
                            pile_lend => $pile_lend,
                            pile_rend => $pile_rend,
                            
                            pile_length => $pile_rend - $pile_lend + 1,

                    
                            new_pile_lend => $pile_lend,
                            new_pile_rend => $pile_rend,
        };
        push (@pile_structs, $pile_struct);
    }
   

    @pile_structs = sort { $a->{pile_lend} <=> $b->{pile_lend} } @pile_structs;

    ## set new pile bounds based on max intron length
    for (my $i = 1; $i <= $#pile_structs; $i++) {
        my $prev_pile_struct = $pile_structs[$i-1];
        my $curr_pile_struct = $pile_structs[$i];

        my $intron_length = $curr_pile_struct->{pile_lend} - $prev_pile_struct->{pile_rend} - 1;
        if ($intron_length > $max_intron_length) {
            $intron_length = $max_intron_length;
        }
        $curr_pile_struct->{new_pile_lend} = $prev_pile_struct->{new_pile_rend} + $intron_length + 1;
        $curr_pile_struct->{new_pile_rend} = $curr_pile_struct->{new_pile_lend} + $curr_pile_struct->{pile_length} - 1;
        
    }

    ## adjust gtf exon coordinates
    
    my $gtf_adj = "";
    my $gene_seq_adj = "";
    
    my $prev_old_pile_rend = 0;
    
    foreach my $pile_struct (@pile_structs) {
        
        my $old_pile_lend = $pile_struct->{pile_lend};
        my $new_pile_lend = $pile_struct->{new_pile_lend};
        my $pile_length = $pile_struct->{pile_length};
        
        my $delta = $old_pile_lend - $new_pile_lend;
        
        ## add intron
        my $intron_len = $old_pile_lend - $prev_old_pile_rend -1;
        my $intron_seq = "";
        if ($prev_old_pile_rend == 0 || $intron_len < $max_intron_length) {
            $intron_seq = substr($gene_seq_region, $prev_old_pile_rend, $intron_len);
        }
        else {
            ## split the difference
            my $left_intron_size = int($max_intron_length/2);
            my $right_intron_size = $max_intron_length - $left_intron_size;
            $left_intron_size -= 5;
            $right_intron_size -= 5; # add 10 Ns at center
            $intron_seq = substr($gene_seq_region, $prev_old_pile_rend, $left_intron_size);
            $intron_seq .= 'N' x 10;
            $intron_seq .= substr($gene_seq_region, $old_pile_lend - 1 - $right_intron_size, $right_intron_size);
            
            if (length($intron_seq) != $max_intron_length) {
                die "Error, intron length is off: " . length($intron_seq) . " vs. $max_intron_length (max)";
            }
        }
        $gene_seq_adj .= $intron_seq;

        foreach my $gtf_row_aref (@{$pile_struct->{pile}}) {
            
            $gtf_row_aref->[3] -= $delta;
            $gtf_row_aref->[4] -= $delta;

            $gtf_adj .= join("\t", @$gtf_row_aref) . "\n";
        }

        my $pile_seq = substr($gene_seq_region, $old_pile_lend -1, $pile_length);
        $gene_seq_adj .= $pile_seq;

        $prev_old_pile_rend = $pile_struct->{pile_rend};
        
    }
    
    ## tack on end of sequence
    $gene_seq_adj .= substr($gene_seq_region, $prev_old_pile_rend);
    
    return($gtf_adj, $gene_seq_adj);
}


    
####
sub set_gtf_scaffold_name {
    my ($scaffold_name, $gtf_text) = @_;

    my $new_gtf = "";
    
    foreach my $line (split(/\n/, $gtf_text)) {
        
        my @x = split(/\t/, $line);
        $x[0] = $scaffold_name;
        
        $x[8] =~ s/transcript_id \"/transcript_id \"$scaffold_name\^/;
        $x[8] =~ s/gene_id \"/gene_id \"$scaffold_name\^/;
        
        $new_gtf .= join("\t", @x) . "\n";
    }
    
    return($new_gtf);
}

####
sub get_gene_contig_gtf {
    my ($gene_gtf, $genome_fasta_file) = @_;

    
    my ($gene_chr, $gene_lend, $gene_rend, $gene_orient, $revised_gene_gtf) = &get_gene_span_info($gene_gtf);
    
    $gene_gtf = $revised_gene_gtf;

    # print STDERR "INFO: gene_chr: $gene_chr, lend: $gene_lend, rend: $gene_rend\n";

    #print STDERR "\n\nGENE_GTF:\n$gene_gtf\n\n";
    
    
    my $seq_region = &get_genomic_region_sequence($genome_fasta_file,
                                                  $gene_chr, 
                                                  $gene_lend - $genome_flank_size, 
                                                  $gene_rend + $genome_flank_size,
                                                  $gene_orient);

    my $gene_contig_gtf = &transform_gtf_coordinates($gene_lend - $genome_flank_size,
                                                     $gene_gtf,
                                                     length($seq_region),
                                                     $gene_orient);
    

    # print STDERR "GENE_contig_gtf:\n$gene_contig_gtf\n\n";

    return($gene_contig_gtf, $seq_region);
    
    
}


#####
sub get_genomic_region_sequence {
    my ($fasta_file, $chr, $lend, $rend, $orient) = @_;

    my $cmd = "samtools faidx $fasta_file $chr:$lend-$rend";
    my $seq = `$cmd`;
    if ($?) {
        die "Error, cmd: $cmd died with ret $?";
    }
    my $header;
    ($header, $seq) = split(/\n/, $seq, 2);
    $seq =~ s/\s//g;
    
    my $seq_len = $rend - $lend + 1;
    if (length($seq) != $seq_len) {
        die "Error, didn't extract required sequence from $fasta_file, $chr, $lend, $rend, instead got seq of length " . length($seq);
    }
    
    if ($orient eq '-') {
        $seq = &reverse_complement($seq);
    }

    return($seq);
}



####
sub extract_gene_gtfs {
    my ($gtf_file, $gene_want_href) = @_;

    my %gene_to_gtf;

    open (my $fh, $gtf_file) or die "Error, cannot open file $gtf_file";
    while (<$fh>) {
        chomp;
        if (/^\#/) { next;}
        unless (/\w/) { next; }
        my $line = $_;
        
        my @x = split(/\t/, $line);
        
        my $chr = $x[0];
        my $feat_type = $x[2];
        my $lend = $x[3];
        my $rend = $x[4];
        my $orient = $x[6];        
        my $info = $x[8];

        unless ($feat_type eq 'exon') { next; } # only exon records of gtf file
        
        my $gene_id = "";
        my $gene_name = "";
        
        # define a gene identifier to use in downstream processes
        # use the ID given by the user
        my $gene_id_use;
        
        if ($info =~ /gene_id \"([^\"]+)\"/) {
            $gene_id = $1;
            $gene_id_use = $gene_id;
        }
        if ($info =~ /gene_name \"([^\"]+)\"/) {
            $gene_name = $1;
            $gene_id_use = $gene_name;
            
            if ($gene_id && $gene_id ne $gene_name) {
                # for viewing purposes
                $info =~ s/$gene_id/$gene_name\^$gene_id/;
                $x[8] = $info;
                $gene_id_use = "$gene_name^$gene_id"; # preferred
            }
        }
        
        
        unless ($gene_want_href->{$gene_id} || $gene_want_href->{$gene_name}) { next; }


        $line = join("\t", @x);
                
        $line .= " FI_gene_label \"$gene_id_use\";";
        
        
        my $orig_info = "$chr,$lend,$rend,$orient";
        $line .= " orig_coord_info \"$orig_info\";\n";
        
        $gene_to_gtf{$gene_id} .= $line;
            
        
        if ($gene_name && $gene_name ne $gene_id) {
            $gene_to_gtf{$gene_name} .= $line;
        }
        
    }
    close $fh;

    return(%gene_to_gtf);
}


####
sub get_gene_span_info {
    my ($gene_gtf_text) = @_;

        
    my ($chr, $min_lend, $max_rend, $orient);

    my $revised_gene_gtf_text = "";

    my @gtf_lines = split(/\n/, $gene_gtf_text);
    foreach my $line (@gtf_lines) {
        my @x = split(/\t/, $line);
        my $scaffold = $x[0];
        my $lend = $x[3];
        my $rend = $x[4];
        ($lend, $rend) = sort {$a<=>$b} ($lend, $rend);
        my $strand = $x[6];
        if (defined $chr) {
            ## check to ensure the rest of the info matches up
            ## if discrepancy, use first found.
            if ($chr ne $scaffold) {
                print STDERR "Error, chr discrepancy in gtf info for $line\n";
                next; 
            }
            if ($orient ne $strand) {
                print STDERR "Error, strand conflict in gtf info for $line\n";
                next;
            }
            if ($lend < $min_lend) {
                $min_lend = $lend;
            }
            if ($rend > $max_rend) {
                $max_rend = $rend;
            }
            
        }
        else {
            ## init
            ($chr, $min_lend, $max_rend, $orient) = ($scaffold, $lend, $rend, $strand);
        }
        $revised_gene_gtf_text .= "$line\n";
    }

    if ($max_rend - $min_lend > 100e6) {
        confess "Error - no gene spans 100M bases in length.... likely problem";
    }
    
    return($chr, $min_lend, $max_rend, $orient, $revised_gene_gtf_text);
}


####
sub transform_gtf_coordinates {
    my ($left_reference_pos, $gene_gtf, $seq_length, $gene_orient) = @_;

    my $new_gtf = "";
    
    foreach my $line (split(/\n/, $gene_gtf)) {
        
        my @fields = split(/\t/, $line);
        my ($lend, $rend) = ($fields[3], $fields[4]);
        
        $lend = $lend - $left_reference_pos + 1;
        $rend = $rend - $left_reference_pos + 1;
        
        if ($gene_orient eq '-') {
            # revcomp the coordinates
            $lend = $seq_length - $lend + 1;
            $rend = $seq_length - $rend + 1;
            ($lend, $rend) = sort {$a<=>$b} ($lend, $rend);
        }

        $fields[3] = $lend;
        $fields[4] = $rend;
        $fields[6] = '+';
        
        $new_gtf .= join("\t", @fields) . "\n";
    }
    
    return($new_gtf);

}

####
sub adjust_gtf_coordinates {
    my ($gtf, $adjustment) = @_;

    my $new_gtf = "";
   
    foreach my $line (split(/\n/, $gtf)) {

        my @fields = split(/\t/, $line);
        $fields[3] += $adjustment;
        $fields[4] += $adjustment;

        $new_gtf .= join("\t", @fields) . "\n";
    }

    return($new_gtf);
}

####
sub process_cmd {
    my ($cmd) = @_;

    print STDERR "CMD: $cmd\n";
    my $ret = system($cmd);
    if ($ret) {
        die "Error, CMD: $cmd died with ret $ret";
    }

    return;
}

