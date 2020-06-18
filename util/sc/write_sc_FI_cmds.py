#!/usr/bin/env python

import sys, os, re
import argparse


def main():


    parser = argparse.ArgumentParser(description="generate FusionInspector commands for each batch of cells", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--batches_list_file", type=str, required=True, help="file containing the list of sample batches")
    parser.add_argument("--genome_lib_dir", type=str, required=True, help="path to the ctat genome lib dir")
    parser.add_argument("--fusions", type=str, required=True, help="targeted fusions file")

    args, unknown_args = parser.parse_known_args()

    batches_list_file = args.batches_list_file
    ctat_genome_lib_dir = args.genome_lib_dir
    fusions_file = args.fusions
    
    FUSION_INSPECTOR_HOME = os.path.dirname( os.path.abspath(__file__)) + "/../.." 
    
    FI_prog = os.path.join(FUSION_INSPECTOR_HOME, "FusionInspector")
    
    with open(batches_list_file, 'rt') as fh:
        for batch in fh:
            batch = batch.rstrip()

            output_dir = batch.replace(".sample_sheet", ".FI.outdir")
            
            cmd = str(FI_prog +
                      " --genome_lib_dir {} ".format(ctat_genome_lib_dir) +
                      " --samples_file {} ".format(batch) +
                      " --fusions {}".format(fusions_file) +
                      " --max_sensitivity " +
                      " --no_FFPM " +
                      " --output_dir {} ".format(output_dir) +
                      " ".join(unknown_args) )

            print(cmd)
    
    sys.exit(0)


if __name__=='__main__':
    main()
