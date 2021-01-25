#!/usr/bin/env python3

import sys, os, re
import csv
from collections import defaultdict
from math import sqrt

MICROH_SIZE = 10
MAX_MICROH_DISTANCE = 10000
MICROH_COUNT_WINDOW_SIZE = 100


def main():

    usage = "usage: {} {} {}\n\n".format(sys.argv[0], "microH.dat.tsv", "fusions.tsv")
    if len(sys.argv) != 3:
        print(usage, file=sys.stderr)
        sys.exit(1)

    microH_dat_filename = sys.argv[1]
    fusions_file = sys.argv[2]
    
    microH_info = defaultdict(list)

    with open(microH_dat_filename, 'rt') as fh:
        for line in fh:
            line = line.rstrip()
            vals = line.split("\t")
            fusion_name = vals[0]
            if vals[1] == "MicroH":
                lend, rend = vals[2], vals[3]
                microH_info[fusion_name].append( [int(lend), int(rend)] )

    
    with open(fusions_file, 'rt') as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fieldnames = list(reader.fieldnames)
        fieldnames.extend(["microh_brkpt_dist", "num_microh_near_brkpt"])


        fusion_name_to_rows = defaultdict(list)
        
        for row in reader:
            fusion_name = row['#FusionName']
            fusion_name_to_rows[fusion_name].append(row)


        writer = csv.DictWriter(sys.stdout, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        
        for fusion_name, preds_list in fusion_name_to_rows.items():
            microH_list = microH_info[fusion_name]
            analyze_fusion_preds(fusion_name, preds_list, microH_list, writer)
            
    
    sys.exit(0)




    
def analyze_fusion_preds(fusion_name, preds_list, microH_list, writer):

    
    for fusion_pred in preds_list:
        left_brkpt = int(fusion_pred['LeftLocalBreakpoint'])
        right_brkpt = int(fusion_pred['RightLocalBreakpoint'])
        
        distances = [MAX_MICROH_DISTANCE] # our local version of infinity
        
        for microH in microH_list:
            distance = get_microh_dist(microH, left_brkpt, right_brkpt)
            distances.append(distance)

        distances = sorted(distances)

        min_dist = distances[0]

        num_within_range = 0
        if min_dist < MAX_MICROH_DISTANCE:
            for dist in distances:
                if dist < MICROH_COUNT_WINDOW_SIZE:
                    num_within_range += 1
                else:
                    break

        fusion_pred['microh_brkpt_dist'] = min_dist
        fusion_pred['num_microh_near_brkpt'] = num_within_range

        writer.writerow(fusion_pred)
        
        
    return

            

def get_microh_dist(microh, geneA_brkpt, geneB_brkpt):

    microh_geneA_lend = microh[0]
    microh_geneA_rend = microh_geneA_lend +  MICROH_SIZE - 1 

    microh_geneB_lend = microh[1]
    microh_geneB_rend = microh_geneB_lend + MICROH_SIZE - 1 
    
    geneA_delta = min( abs(microh_geneA_lend - geneA_brkpt), abs(microh_geneA_rend - geneA_brkpt) )
    geneB_delta = min( abs(microh_geneB_lend - geneB_brkpt), abs(microh_geneB_rend - geneB_brkpt) )

    distance = int( sqrt( geneA_delta**2 + geneB_delta**2) )
    
    return distance


    

    

if __name__=='__main__':
    main()
