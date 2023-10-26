#!/usr/bin/env python3
# encoding: utf-8


import argparse
import csv
import json
import gzip
import os
import sys
import logging

logging.basicConfig(stream=sys.stderr, level=logging.INFO)
logger = logging.getLogger(__name__)


# Change headers to more readable headers for the data table
convert_header_to_eng = {
    "#FusionName" : "Fusion",
    "num_LR" : "# Long Reads",
    "JunctionReadCount" : "Junction Reads",
    "SpanningFragCount" : "Spanning Fragments",
    "FFPM" : "Expr Level (FFPM)",
    "SpliceType" : "Splice Type",
    "LeftGene" : "Left Gene",
    "RightGene" : "Right Gene",
    "LeftBreakpoint" : "Left Breakpoint",
    "RightBreakpoint" : "Right Breakpoint",
    "annots" : "Annotations"
}


arguments = argparse.ArgumentParser( prog = "Fusion Inspector JSON Maker",
                                     description = "Makes a JSON file for a directory of results from fusion inspector",
                                     formatter_class = argparse.ArgumentDefaultsHelpFormatter )

arguments.add_argument( "--fusion_inspector_directory",
                        dest="fusion_inspector_directory",
                        required=True,
                        type=str,
                        help = "The input directory to create the json from; this folder should " +
                        "contain the finspector.fusion_predictions.final.abridged file")

arguments.add_argument("--json_outfile",
                       dest="output_json_file",
                       required=True,
                       type=str,
                       help = "The output json file to create" )


arguments.add_argument("--roi_outfile",
                       required=True,
                       type=str,
                       help="bed file indicating fusion breakpoints as roi")
                       


arguments.add_argument("--file_prefix",
                       dest = "file_prefix",
                       required=True,
                       type=str,
                       help="prefix to FusionInspector output files")


args = arguments.parse_args()

# Make sure the input directory is an absolute path
absolute_fusion_directory = os.path.abspath( args.fusion_inspector_directory )

file_prefix = args.file_prefix

# Include Trinity related files
C_STR_INCLUDE_TRINITY_BED = file_prefix + ".gmap_trinity_GG.fusions.gff3.bed.sorted.bed"
C_STR_INCLUDE_TRINITY_BED_GZ = file_prefix + ".gmap_trinity_GG.fusions.gff3.bed.sorted.bed.gz"


# Dict to be translated to JSON object
dict_json = {"fusions" : [] }


fusions_table = os.path.join( absolute_fusion_directory, args.file_prefix + ".FusionInspector.fusions.abridged.tsv")


fusion_breakpoint_rois = list()

# Make fusion detail
with open(fusions_table, "rt" ) as fh:
    fusion_detail_reader = csv.DictReader(fh, delimiter="\t")

    column_headers = list(fusion_detail_reader.fieldnames)
    for column_header in column_headers:
        if column_header not in convert_header_to_eng:
            logger.warn(f"Missing text conversion for header field: {column_header}")

    # Parse fusion annotation information
    for row in fusion_detail_reader:
        fusion_detail_current = {}
        
        for column_name in column_headers:
            if column_name in convert_header_to_eng:
                fusion_detail_current[ convert_header_to_eng[column_name] ] = row[column_name]

                        
        dict_json["fusions"].append( fusion_detail_current )


        # get roi breakpoint info:
        fusion_breakpoint_rois.append([row['#FusionName'],
                                       row['LeftLocalBreakpoint'],
                                       str(int(row['RightLocalBreakpoint']) -1)])

# Store as a json object
with open( args.output_json_file, "w" ) as write_json:
    write_json.write( json.dumps( dict_json, sort_keys=True, indent= 2 ) )


# write roi bed file
with open(args.roi_outfile, "wt") as ofh:
    for roi in fusion_breakpoint_rois:
        print("\t".join(roi), file=ofh)
        


sys.exit(0)

