#!/usr/bin/env python3

import os, re, sys
import argparse
import subprocess
import csv

import logging

FORMAT = "%(asctime)-15s %(levelname)s %(module)s.%(name)s.%(funcName)s at %(lineno)d :\n\t%(message)s\n"
global logger
logger = logging.getLogger()
logging.basicConfig(level=logging.DEBUG)


UTILDIR = os.path.dirname(__file__)
MISCDIR = os.path.join(UTILDIR, "misc")


RG_OBJ_FILE = os.path.join(MISCDIR, "data", "ranger.rg_obj.rds")


def main():


    arg_parser = argparse.ArgumentParser()

    arg_parser.add_argument("--FI_results_tsv",
                            type=str,
                            required=True,
                            help="FI fusion results")

    arg_parser.add_argument(
        "--genome_lib_dir",
        dest="genome_lib_dir",
        type=str,
        default=os.environ.get("CTAT_GENOME_LIB"),
        help="genome lib directory - see http://FusionFilter.github.io for details.  Uses env var CTAT_GENOME_LIB as default",
        )
    
    
    args_parsed = arg_parser.parse_args()

    FI_results_file = args_parsed.FI_results_tsv
        
    genome_lib_dir = args_parsed.genome_lib_dir
    if genome_lib_dir is None:
        raise RuntimeError(
        "Error, must specify --genome_lib_dir or set env var CTAT_GENOME_LIB"
            )

    genome_lib_dir = os.path.abspath(genome_lib_dir)

    ## begin work:

    
    # prep for cluster prediction
    cmdstr = str(" ".join([ os.path.join(MISCDIR, "prep_data_for_cosmic-like_pred.py"),
                            FI_results_file,
                            "> {}.prepped".format(FI_results_file) ]) )
    run_cmd(cmdstr)


    # predict cosmic-like cluster
    cmdstr = str(" ".join([ os.path.join(MISCDIR, "predict_cosmic_like_fusion_cluster.R"),
                            "--fusions {}.prepped".format(FI_results_file),
                            "--ranger {}".format(RG_OBJ_FILE),
                            "--output {}.wPreds.tsv".format(FI_results_file)]) )

    run_cmd(cmdstr)

    logger.info("Done")
    

    sys.exit(0)



def run_cmd(cmdstr):

    logger.info("CMD: {}".format(cmdstr))

    subprocess.check_call(cmdstr, shell=True)

    return



def write_fusions_list_file(FI_results_file, fusions_list_file):

    fusions = set()

    with open(FI_results_file) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            fusion_name = row['#FusionName']
            fusions.add(fusion_name)

    with open(fusions_list_file, "wt") as ofh:
        for fusion in fusions:
            print(fusion, file=ofh)

    return





if __name__=='__main__':
    main()




