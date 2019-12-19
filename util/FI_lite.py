#!/usr/bin/env python
# encoding: utf-8


import os, re, sys
import argparse
import subprocess

sys.path.insert(0, os.path.sep.join([os.path.dirname(os.path.realpath(__file__)), "../PyLib"]))
from Pipeliner import Pipeliner, Command

import logging
FORMAT = "%(asctime)-15s %(levelname)s %(module)s.%(name)s.%(funcName)s at %(lineno)d :\n\t%(message)s\n"
global logger
logger = logging.getLogger()
logging.basicConfig(filename='FusionInspector.log', format=FORMAT, filemode='w', level=logging.DEBUG)
# add a new Handler to print all INFO and above messages to stdout
ch = logging.StreamHandler(sys.stdout)
ch.setLevel(logging.INFO)
logger.addHandler(ch)

UTILDIR = os.path.dirname(os.path.realpath(__file__))

                          
def main():

    arg_parser = argparse.ArgumentParser(description="FI lite - for debugging purposes",
                                        formatter_class=argparse.RawTextHelpFormatter)
    
    arg_parser.add_argument("--fi_gtf", type=str, required=True, help="finspector gtf file")
    arg_parser.add_argument("--bam", type=str, required=True, help="alignments to fi contigs bam file")
    arg_parser.add_argument("--min_per_id", dest="min_per_id", type=int, required=False, default=96,
                            help='minimum percent identity for a fusion-supporting read alignment')
    
    args = arg_parser.parse_args()

    args_parsed = args
    bam_file = args.bam
    mergedContig_gtf_filename = args.fi_gtf

    chckpts_dir = "__chckpts.{}".format(os.getpid())
    if not os.path.exists(chckpts_dir):
        os.makedirs(chckpts_dir)
    
    pipeliner = Pipeliner(chckpts_dir)
    
    ## extract the fusion JUNCTION reads
    fusion_junction_reads_sam_file = bam_file + ".fusion_junc_reads.sam"
    fusion_junction_info_file = bam_file + ".fusion_junction_info"
        
    cmdstr = str(os.sep.join([UTILDIR, "get_fusion_JUNCTION_reads_from_fusion_contig_bam.pl"]) +
                 " --gtf_file " + mergedContig_gtf_filename +
                 " --MIN_ALIGN_PER_ID " + str(args_parsed.min_per_id) +
                 " --bam " + bam_file + " > " + fusion_junction_reads_sam_file)
    
    pipeliner.add_commands([Command(cmdstr, "get_fusion_JUNCTION_reads_from_bam.ok")]) 

    
    ## extract the fusion SPANNING reads
    fusion_spanning_reads_sam_file = bam_file + ".fusion_span_reads.sam"
    fusion_spanning_reads_info_file = bam_file + ".fusion_spanning_info"
    
    cmdstr = str(os.sep.join([UTILDIR, "get_fusion_SPANNING_reads_from_bam.from_chim_summary.pl"]) + 
                 " --gtf_file " + mergedContig_gtf_filename +
                 " --MIN_ALIGN_PER_ID " + str(args_parsed.min_per_id) +
                 " --bam " + bam_file +
                 " --junction_info " + fusion_junction_info_file + 
                 " > " + fusion_spanning_reads_sam_file)
    
    pipeliner.add_commands([Command(cmdstr, "get_fusion_SPANNING_reads_from_bam.ok")]) 

    
    ## coalesce junction and spanning info
    fusion_summary_file = "FI_lite.fusions.unfiltered.tsv"
    cmdstr = str( os.sep.join([UTILDIR, "coalesce_junction_and_spanning_info.pl"]) + " " +
                 fusion_junction_info_file + " " +
                 fusion_spanning_reads_info_file + " " +
                      " 1 " + # FAR pseudocount
                      " > " + fusion_summary_file )
    pipeliner.add_commands([Command(cmdstr, "coalesce_junc_n_span.ok")]) 


    abridged_summary_file = "FI_lite.fusions.unfiltered.abridged.tsv"
    cmdstr = str(UTILDIR + "/column_exclusions.pl " + fusion_summary_file +
                 " JunctionReads,SpanningFrags,CounterFusionLeftReads,CounterFusionRightReads " +
                 " > " + abridged_summary_file );
    pipeliner.add_commands([Command(cmdstr, "final.abridged.ok")]) 
    
    pipeliner.run()

    sys.exit(0)

if __name__ == '__main__':
    main()
