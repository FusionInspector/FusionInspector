#!/usr/bin/env python3
# encoding: utf-8

import os, sys, re
import logging
import argparse
import pysam

logging.basicConfig(level=logging.INFO, 
                    format='%(asctime)s : %(levelname)s : %(message)s',
                    datefmt='%H:%M:%S')
logger = logging.getLogger(__name__)

def main():

    parser = argparse.ArgumentParser(description="mark duplicates in bam", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--input_bam", "-i", dest="input_bam", required=True, type=str,
                        help="input bam file, coordinate sorted")

    parser.add_argument("--output_bam", "-o", dest="output_bam", required=True, type=str,
                        help="output bam file")

    parser.add_argument("--remove_dups", "-r", dest='remove_dups', action='store_true', default=False,
                        help='instead of marking duplicates, just remove them')
    
    args = parser.parse_args()

    input_bam_filename = args.input_bam
    output_bam_filename = args.output_bam
    remove_dups_flag = args.remove_dups
    
    bamreader = pysam.AlignmentFile(input_bam_filename, "rb")

    if ( ( (not 'SO' in bamreader.header.as_dict()['HD']) )
        or
        bamreader.header.as_dict()['HD']['SO'] != 'coordinate') :
        raise RuntimeError("Error, file: {} must be coordinate sorted".format(input_bam_filename))
    
    bamwriter =  pysam.AlignmentFile(output_bam_filename, "wb", template=bamreader)
    

    # KISS: just use the read and mate starting points
    
    prev_start = -1
    prev_chrom = None
    queued_duplicate_reads = dict()
    current_pos_mate_coords = set()

    def reinit_current_contig():
        nonlocal current_pos_mate_coords
        current_pos_mate_coords.clear()

    def reinit_new_contig():
        reinit_current_contig()
        nonlocal queued_duplicate_reads
        queued_duplicate_reads.clear()
        prev_start = -1
        prev_chrom = None
        
        
    for read in bamreader.fetch():
        chrom = bamreader.get_reference_name(read.reference_id)
        start = read.reference_start
        read_name = read.query_name

        if read.is_secondary:
            continue
        
        if chrom != prev_chrom:
            reinit_new_contig()
        elif start != prev_start:
            reinit_current_contig()
        
        mate_start = read.next_reference_start

        duplicate_flag = False
        
        if read_name in queued_duplicate_reads and queued_duplicate_reads[read_name] == start:
            # mark as duplicate
            duplicate_flag = True
            del queued_duplicate_reads[read_name]
        
        elif mate_start > start and mate_start in current_pos_mate_coords:
            # mark this read as a duplicate
            queued_duplicate_reads[read_name] = mate_start
            duplicate_flag = True

        else:
            # not a duplicate
            # store mate coord in case we find others w/ similar read and mate starts
            current_pos_mate_coords.add(mate_start)


        # output read alignment.
        if duplicate_flag:
            read.is_duplicate = True

        if (not duplicate_flag) or (duplicate_flag and not remove_dups_flag):
            bamwriter.write(read)
            
        prev_chrom = chrom
        prev_start = start


    logger.info("Done.")
    
    sys.exit(0)


if __name__=='__main__':
    main()
