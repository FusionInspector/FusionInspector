#!/usr/bin/env python

import argparse
import csv
import json
import gzip
import os

# Argument
C_ARG_INCLUDE_TRINITY = "--include_Trinity"

# Include Trinity related files
C_STR_INCLUDE_TRINITY_BED = "finspector.gmap_trinity_GG.fusions.gff3.bed.sorted.bed"
C_STR_INCLUDE_TRINITY_BED_GZ = "finspector.gmap_trinity_GG.fusions.gff3.bed.sorted.bed.gz"

# Constants keys in JSON file
C_REFERENCE = "reference"
C_REFERENCE_INDEX = "referenceIndex"
C_CYTOBAND = "cytoband"
C_REFERENCE_BED = "referenceBed"
C_JUNCTION_SPANNING = "junctionSpanning"
C_TRINITY_BED = "trinityBed"
C_JUNCTION_READS = "junctionReads"
C_JUNCTION_READS_BAM = "junctionReadsBam"
C_JUNCTION_READS_BAI = "junctionReadsBai"
C_SPANNING_READS = "spanningReads"
C_SPANNING_READS_BAM = "spanningReadsBam"
C_SPANNING_READS_BAI = "spanningReadsBai"
C_SAMPLE_NAME = "sampleName"
C_FUSION_DETAIL = "fusions"
C_FUSION_DETAIL_NAME = "Name"
C_FUSION_DETAIL_LEFT_CHR = "Left Chr"
C_FUSION_DETAIL_LEFT_POS = "Left Pos"
C_FUSION_DETAIL_LEFT_STRAND = "Left Strand"
C_FUSION_DETAIL_RIGHT_CHR = "Right Chr"
C_FUSION_DETAIL_RIGHT_POS = "Right Pos"
C_FUSION_DETAIL_RIGHT_STRAND = "Right Strand"
C_FUSION_DETAIL_LEFT_GENE = "Left Gene"
C_FUSION_DETAIL_RIGHT_GENE = "Right Gene"
C_FUSION_DETAIL_SPLICE_TYPE = "Splice Type"
C_FUSION_DETAIL_JUNCTION_READS = "Junction Reads"
C_FUSION_DETAIL_SPANNING_FRAGMENTS = "Spanning Frags"
C_FUSION_INSPECTOR_DIR = "FUSION_INSPECTOR"

# Constants from Fusion Inspector (in output file that is parsed into the fusion detail)
C_PRED_FUSION_NAME = "#fusion_name"
C_PRED_FUSION_NAME_ENG = "Fusion"
C_PRED_JUNCTION_READS = "JunctionReads"
C_PRED_JUNCTION_READS_ENG = "Junction Reads"
C_PRED_SPANNING_FRAGS = "SpanningFrags"
C_PRED_SPANNING_FRAGS_ENG = "Spanning Fragments"
C_PRED_SPLICE_TYPE = "Splice_type"
C_PRED_SPLICE_TYPE_ENG = "Splice Type"
C_PRED_LEFT_GENE = "LeftGene"
C_PRED_LEFT_GENE_ENG = "Left Gene"
C_PRED_RIGHT_GENE = "RightGene"
C_PRED_RIGHT_GENE_ENG = "Right Gene"
C_PRED_LEFT_BREAK_POINT = "LeftBreakpoint"
C_PRED_LEFT_BREAK_POINT_ENG = "Left Breakpoint"
C_PRED_RIGHT_BREAK_POINT = "RightBreakpoint"
C_PRED_RIGHT_BREAK_POINT_ENG = "Right Breakpoint"

# Convert splice type token to a more readable form
C_SPLICE_REFERENCE = "INCL_NON_REF_SPLICE"
C_SPLICE_REFERENCE_ENG = "Includes Reference"
C_SPLICE_NOT_REFERENCE = "ONLY_REF_SPLICE"
C_SPLICE_NOT_REFERENCE_ENG = "DOES NOT Include Reference"

# Change splice info to plain english 
convert_splice_type_to_eng = {
    C_SPLICE_REFERENCE : C_SPLICE_REFERENCE_ENG,
    C_SPLICE_NOT_REFERENCE : C_SPLICE_NOT_REFERENCE_ENG
}

# Change headers to more readable headers for the data table
convert_header_to_eng = {
    C_PRED_FUSION_NAME : C_PRED_FUSION_NAME_ENG,
    C_PRED_JUNCTION_READS : C_PRED_JUNCTION_READS_ENG,
    C_PRED_SPANNING_FRAGS : C_PRED_SPANNING_FRAGS_ENG,
    C_PRED_SPLICE_TYPE : C_PRED_SPLICE_TYPE_ENG,
    C_PRED_LEFT_GENE : C_PRED_LEFT_GENE_ENG,
    C_PRED_RIGHT_GENE : C_PRED_RIGHT_GENE_ENG,
    C_PRED_LEFT_BREAK_POINT : C_PRED_LEFT_BREAK_POINT_ENG,
    C_PRED_RIGHT_BREAK_POINT : C_PRED_RIGHT_BREAK_POINT_ENG
}

arguments = argparse.ArgumentParser( prog = "Fusion Inspector JSON Maker", description = "Makes a JSON file for a directory of results from fusion inspector", formatter_class = argparse.ArgumentDefaultsHelpFormatter )
arguments.add_argument( C_ARG_INCLUDE_TRINITY, dest="f_include_trinity", action="store_true", help = "Adds the gmap_trinity_GG bed file to the json for viewing." )
arguments.add_argument( dest="fusion_inspector_directory", help = "The input directory to create the json from; this folder should contain the " + C_FUSION_INSPECTOR_DIR + " directory." )
arguments.add_argument( dest="output_json_file", help = "The output json file to create" )
args = arguments.parse_args()

# Make sure the input directory is an absolute path
absolute_fusion_directory = os.path.abspath( args.fusion_inspector_directory )

# Dict to be translated to JSON object
dict_json = {}
dict_json[ C_REFERENCE ] = "finspector.fa"
dict_json[ C_REFERENCE_INDEX ] = "finspector.fa.fai"
dict_json[ C_CYTOBAND ] = "cytoBand.txt"
dict_json[ C_REFERENCE_BED ] = "finspector.bed"
dict_json[ C_JUNCTION_SPANNING ] = "finspector.igv.FusionJuncSpan"
dict_json[ C_TRINITY_BED ] = C_STR_INCLUDE_TRINITY_BED if args.f_include_trinity else "NA"
dict_json[ C_JUNCTION_READS ] = "finspector.junction_reads.bam.bed.sorted.bed.gz"
dict_json[ C_JUNCTION_READS_BAM ] = "finspector.junction_reads.bam"
dict_json[ C_JUNCTION_READS_BAI ] = "finspector.junction_reads.bam.bai"
dict_json[ C_SPANNING_READS ] = "finspector.spanning_reads.bam.bed.sorted.bed.gz"
dict_json[ C_SPANNING_READS_BAM ] = "finspector.spanning_reads.bam"
dict_json[ C_SPANNING_READS_BAI ] = "finspector.spanning_reads.bam.bai"
dict_json[ C_SAMPLE_NAME ] = os.path.basename( args.fusion_inspector_directory ) 
dict_json[ C_FUSION_DETAIL ] = []

# Uncompress bed file for Galaxy (which does not like compressed files).
if args.f_include_trinity:
    str_trinity_bed_file = os.path.join( absolute_fusion_directory, C_FUSION_INSPECTOR_DIR, C_STR_INCLUDE_TRINITY_BED )
    str_trinity_bed_file_compressed = os.path.join( absolute_fusion_directory, C_FUSION_INSPECTOR_DIR, C_STR_INCLUDE_TRINITY_BED_GZ )
    if not os.path.exists( str_trinity_bed_file ):
        if os.path.exists( str_trinity_bed_file_compressed ):
            with gzip.open( str_trinity_bed_file_compressed ) as hndl_compressed_file:
                with open( str_trinity_bed_file, "w" ) as hndl_uncompressed_file:
                    for str_compressed_line in hndl_compressed_file:
                        hndl_uncompressed_file.write( str_compressed_line ) 
        else:
            print C_ARG_INCLUDE_TRINITY + " was given but one of the following files are expected to exist and did not:" + " or ".join([ C_STR_INCLUDE_TRINITY_BED, C_STR_INCLUDE_TRINITY_BED_GZ ])
            exit( -1 )

# Make fusion detail
with open( os.path.join( absolute_fusion_directory, C_FUSION_INSPECTOR_DIR, "finspector.fusion_predictions.final.abridged" ), "r" ) as fusion_detail_contents:
    fusion_detail_parsed = csv.reader( fusion_detail_contents, delimiter = "\t" )

    # Get indices of the columns of interest, if they do not exist error.
    detail_key_to_index = {}
    header_line = fusion_detail_parsed.next()
    header_keys = [ C_PRED_FUSION_NAME, C_PRED_JUNCTION_READS, C_PRED_SPANNING_FRAGS, C_PRED_SPLICE_TYPE,
                      C_PRED_LEFT_GENE, C_PRED_RIGHT_GENE, C_PRED_LEFT_BREAK_POINT, C_PRED_RIGHT_BREAK_POINT ]

    # Get indices to shuffle incoming header to the given order in header_keys
    for header_key in header_keys:
        detail_key_to_index[ header_key] = header_line.index( header_key )
        if( detail_key_to_index[ header_key] < 0 ):
            print( "Could not find the index for " + header_key + "." )

    # Parse fusion annotation information
    for fusion_detail_line in fusion_detail_parsed:
        fusion_detail_current = {}

        # Parse compound annotations
        for header_key in header_keys:
            if header_key == C_PRED_LEFT_BREAK_POINT:
                break_point_entries = fusion_detail_line[ detail_key_to_index[ header_key ] ].split(":")
                fusion_detail_current[ C_FUSION_DETAIL_LEFT_CHR ] = break_point_entries[ 0 ]
                fusion_detail_current[ C_FUSION_DETAIL_LEFT_POS ] = break_point_entries[ 1 ]
                fusion_detail_current[ C_FUSION_DETAIL_LEFT_STRAND ] = break_point_entries[ 2 ]
            elif header_key == C_PRED_RIGHT_BREAK_POINT:
                break_point_entries = fusion_detail_line[ detail_key_to_index[ header_key ] ].split(":")
                fusion_detail_current[ C_FUSION_DETAIL_RIGHT_CHR ] = break_point_entries[ 0 ]
                fusion_detail_current[ C_FUSION_DETAIL_RIGHT_POS ] = break_point_entries[ 1 ]
                fusion_detail_current[ C_FUSION_DETAIL_RIGHT_STRAND ] = break_point_entries[ 2 ]
            elif header_key == C_PRED_SPLICE_TYPE:
                fusion_detail_current[ convert_header_to_eng[ header_key ] ] = convert_splice_type_to_eng[ fusion_detail_line[ detail_key_to_index[ header_key ] ] ]
            elif header_key == C_PRED_FUSION_NAME:
                fusion_detail_current[ convert_header_to_eng[ header_key ] ] = fusion_detail_line[ detail_key_to_index[ header_key ] ]
            else:
                fusion_detail_current[ convert_header_to_eng[ header_key ] ] = fusion_detail_line[ detail_key_to_index[ header_key ] ]

    dict_json[ C_FUSION_DETAIL ].append( fusion_detail_current )

# Store as a json object
with open( args.output_json_file, "w" ) as write_json:
    write_json.write( json.dumps( dict_json, sort_keys=True, indent= 2 ) )
