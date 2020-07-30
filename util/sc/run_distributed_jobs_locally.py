#!/usr/bin/env python

import sys, os, re
import argparse
import subprocess
import logging

logging.basicConfig(stream=sys.stderr, level=logging.INFO)
logger = logging.getLogger(__name__)



def main():

    parser = argparse.ArgumentParser(
        description="generate FusionInspector commands for each batch of cells",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("--cmds_file", type=str, required=True, help="file containing the list of FusionInspector commands")
    parser.add_argument("--num_parallel_exec", type=int, required=True, help="number of FusionInspector commands to run in parallel")
    parser.add_argument("--genome_lib_dir", type=str, required=True, help="path to ctat genome lib dir")


    args = parser.parse_args()

    cmds_file = args.cmds_file
    num_parallel = args.num_parallel_exec
    genome_lib_dir = args.genome_lib_dir


    ensure_parafly_installed()
    

    cmd = "ParaFly -c {} -CPU {} -vv -max_retry 1".format(cmds_file, num_parallel)
    subprocess.check_call(cmd, shell=True)

    logger.info("-done running parallel commands.")


    logger.info("done")
    
    sys.exit(0)
    

def ensure_parafly_installed():

    try:
        subprocess.check_call("which ParaFly", shell=True)
    except:
        raise RuntimeError("Erorr, cannot find ParaFly installed and available in the PATH setting.  Be sure to install ParaFly before running")

    return


    
if __name__=="__main__":
    main()

