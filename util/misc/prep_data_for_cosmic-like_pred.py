#!/usr/bin/env python3

import sys, os, re
import csv



def main():


    usage = "usage: {} fusions.wMicroH.tsv\n\n".format(sys.argv[0])
    if len(sys.argv) < 2:
        print(usage, file=sys.stderr)
        sys.exit(1)


    fusions_file = sys.argv[1]


    with open(fusions_file, 'rt') as fh:
        reader = csv.DictReader(fh, delimiter="\t")
    
        colnames = list(reader.fieldnames)

        colnames.extend(["annot_splice", "consensus_splice", "left_counter_ffpm", "right_counter_ffpm"])
            
        writer = csv.DictWriter(sys.stdout, fieldnames=colnames, delimiter="\t")
    
        writer.writeheader()
    

        for row in reader:
            splice_type = 1 if (row['SpliceType'] == "ONLY_REF_SPLICE") else 0
            ffpm = float(row['FFPM'])
            left_dinuc = row['LeftBreakDinuc'].upper()
            right_dinuc = row['RightBreakDinuc'].upper()

            dinuc_pair = left_dinuc + ":" + right_dinuc
                
            consensus_splice = 1 if dinuc_pair in ["GT:AG", "GC:AG"] else 0
            
            num_tot_reads = float(row['est_J']) + float(row['est_S'])
            
            ffpm_scaling_factor = ffpm /num_tot_reads  

            left_counter_ffpm = int(row['NumCounterFusionLeft']) * ffpm_scaling_factor
            right_counter_ffpm = int(row['NumCounterFusionRight']) * ffpm_scaling_factor
            
            row["annot_splice"] = splice_type
            row["consensus_splice"] = consensus_splice
            row["left_counter_ffpm"] = left_counter_ffpm
            row["right_counter_ffpm"] = right_counter_ffpm
            
            
            writer.writerow(row)
            
            

    sys.exit(0)


if __name__=='__main__':
    main()
