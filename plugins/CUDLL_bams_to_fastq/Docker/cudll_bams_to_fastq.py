#!/usr/bin/env python3

"""
Extract primary reads from two CUDLL BAM files and merge into a single FASTQ.
Reads from the supplemental BAM take priority; reads with the same name in the
main BAM are excluded.
"""

import argparse
import sys
import gzip
import pysam
from collections import defaultdict


def extract_read_names(bam_path):
    """Extract all primary read names from a BAM file."""
    read_names = set()
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam:
            # Only process primary alignments
            if not read.is_secondary and not read.is_supplementary:
                read_names.add(read.query_name)
    return read_names


def write_fastq_record(outfile, read):
    """Write a single read to FASTQ format."""
    # Get the sequence and quality
    seq = read.query_sequence
    qual = pysam.qualities_to_qualitystring(read.query_qualities)

    # Handle reverse complement if needed
    if read.is_reverse:
        # pysam returns the sequence in forward orientation by default
        # query_sequence is already oriented as it appears in the FASTQ
        pass

    # Determine read pair suffix
    if read.is_paired:
        if read.is_read1:
            suffix = "/1"
        elif read.is_read2:
            suffix = "/2"
        else:
            suffix = ""
    else:
        suffix = ""

    # Write FASTQ record
    outfile.write(f"@{read.query_name}{suffix}\n")
    outfile.write(f"{seq}\n")
    outfile.write("+\n")
    outfile.write(f"{qual}\n")


def process_bams(main_bam, supp_bam, output_fastq):
    """
    Process two BAM files and create a merged FASTQ.
    Priority is given to reads from supp_bam.
    """
    print(f"Extracting read names from supplemental BAM: {supp_bam}", file=sys.stderr)
    supp_read_names = extract_read_names(supp_bam)
    print(f"Found {len(supp_read_names)} unique read names in supplemental BAM", file=sys.stderr)

    # Open output file
    with gzip.open(output_fastq, "wt") as outfile:
        # First, write all primary reads from supplemental BAM
        print(f"Writing reads from supplemental BAM...", file=sys.stderr)
        supp_count = 0
        with pysam.AlignmentFile(supp_bam, "rb") as bam:
            for read in bam:
                if not read.is_secondary and not read.is_supplementary:
                    write_fastq_record(outfile, read)
                    supp_count += 1

        print(f"Wrote {supp_count} reads from supplemental BAM", file=sys.stderr)

        # Then, write primary reads from main BAM that are NOT in supplemental BAM
        print(f"Writing non-overlapping reads from main BAM: {main_bam}", file=sys.stderr)
        main_count = 0
        excluded_count = 0
        with pysam.AlignmentFile(main_bam, "rb") as bam:
            for read in bam:
                if not read.is_secondary and not read.is_supplementary:
                    if read.query_name not in supp_read_names:
                        write_fastq_record(outfile, read)
                        main_count += 1
                    else:
                        excluded_count += 1

        print(f"Wrote {main_count} reads from main BAM", file=sys.stderr)
        print(f"Excluded {excluded_count} reads from main BAM (present in supplemental)", file=sys.stderr)
        print(f"Total reads written: {supp_count + main_count}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description="Extract primary reads from two CUDLL BAM files and merge into FASTQ"
    )
    parser.add_argument(
        "--main-bam",
        required=True,
        help="Main CUDLL BAM file"
    )
    parser.add_argument(
        "--supp-bam",
        required=True,
        help="Supplemental CUDLL BAM file (takes priority)"
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output FASTQ.gz file"
    )

    args = parser.parse_args()

    print(f"Processing CUDLL BAM files...", file=sys.stderr)
    print(f"  Main BAM: {args.main_bam}", file=sys.stderr)
    print(f"  Supp BAM: {args.supp_bam}", file=sys.stderr)
    print(f"  Output: {args.output}", file=sys.stderr)

    process_bams(args.main_bam, args.supp_bam, args.output)

    print("Done!", file=sys.stderr)


if __name__ == "__main__":
    main()
