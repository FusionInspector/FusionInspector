#!/usr/bin/env python3

"""
Extract primary reads from one or two CUDLL BAM files and create a table of read
names with their cell barcodes (CB) and UMIs (XM).
When provided, reads from the supplemental BAM take priority; reads with the
same name in the main BAM are excluded.
"""

import argparse
import sys
import gzip
import pysam
import time


def extract_read_names(bam_path):
    """Extract all primary read names from a BAM file."""
    read_names = set()
    start_time = time.time()
    total_reads = 0
    report_interval = 1000000  # Report every 1M reads

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam:
            total_reads += 1
            # Only process primary alignments
            if not read.is_secondary and not read.is_supplementary:
                read_names.add(read.query_name)

            # Progress reporting
            if total_reads % report_interval == 0:
                elapsed = time.time() - start_time
                rate = total_reads / elapsed
                print(f"  Progress: {total_reads:,} reads processed ({rate:,.0f} reads/sec)",
                      file=sys.stderr, flush=True)

    elapsed = time.time() - start_time
    print(f"  Completed: {total_reads:,} total reads in {elapsed:.1f}s",
          file=sys.stderr, flush=True)
    return read_names


def write_cb_umi_record(outfile, read, cb_tag, umi_tag):
    """Write a single read's CB/UMI info to the output file."""
    read_name = read.query_name

    # Extract cell barcode and UMI from tags
    try:
        cell_barcode = read.get_tag(cb_tag)
    except KeyError:
        cell_barcode = "NA"

    try:
        umi = read.get_tag(umi_tag)
    except KeyError:
        umi = "NA"

    # Write tab-delimited record
    outfile.write(f"{read_name}\t{cell_barcode}\t{umi}\n")


def process_bams(main_bam, supp_bam, output_file, cb_tag, umi_tag):
    """
    Process one or two BAM files and create a CB/UMI table.
    Priority is given to reads from supp_bam when provided.
    """
    supp_read_names = set()
    supp_count = 0

    if supp_bam:
        print(f"Extracting read names from supplemental BAM: {supp_bam}", file=sys.stderr)
        supp_read_names = extract_read_names(supp_bam)
        print(f"Found {len(supp_read_names)} unique read names in supplemental BAM", file=sys.stderr)
    else:
        print("No supplemental BAM provided; processing only the main BAM", file=sys.stderr)

    # Open output file
    opener = gzip.open if output_file.endswith('.gz') else open
    with opener(output_file, "wt") as outfile:
        # Write header
        outfile.write(f"read_name\tcell_barcode\tUMI\n")

        if supp_bam:
            # First, write all primary reads from supplemental BAM
            print("Writing reads from supplemental BAM...", file=sys.stderr)
            supp_missing_cb = 0
            supp_missing_umi = 0
            start_time = time.time()
            report_interval = 1000000  # Report every 1M reads

            with pysam.AlignmentFile(supp_bam, "rb") as bam:
                for read in bam:
                    if not read.is_secondary and not read.is_supplementary:
                        write_cb_umi_record(outfile, read, cb_tag, umi_tag)
                        supp_count += 1
                        # Track missing tags
                        if not read.has_tag(cb_tag):
                            supp_missing_cb += 1
                        if not read.has_tag(umi_tag):
                            supp_missing_umi += 1

                        # Progress reporting
                        if supp_count % report_interval == 0:
                            elapsed = time.time() - start_time
                            rate = supp_count / elapsed
                            print(f"  Progress: {supp_count:,} reads written ({rate:,.0f} reads/sec)",
                                  file=sys.stderr, flush=True)

            elapsed = time.time() - start_time
            print(f"Wrote {supp_count:,} reads from supplemental BAM in {elapsed:.1f}s", file=sys.stderr)
            if supp_missing_cb > 0:
                print(f"  Warning: {supp_missing_cb} reads missing {cb_tag} tag", file=sys.stderr)
            if supp_missing_umi > 0:
                print(f"  Warning: {supp_missing_umi} reads missing {umi_tag} tag", file=sys.stderr)

        # Then, write primary reads from main BAM that are NOT in supplemental BAM
        print(f"Writing non-overlapping reads from main BAM: {main_bam}", file=sys.stderr)
        main_count = 0
        excluded_count = 0
        main_missing_cb = 0
        main_missing_umi = 0
        main_total_processed = 0
        start_time = time.time()

        with pysam.AlignmentFile(main_bam, "rb") as bam:
            for read in bam:
                if not read.is_secondary and not read.is_supplementary:
                    main_total_processed += 1
                    if read.query_name not in supp_read_names:
                        write_cb_umi_record(outfile, read, cb_tag, umi_tag)
                        main_count += 1
                        # Track missing tags
                        if not read.has_tag(cb_tag):
                            main_missing_cb += 1
                        if not read.has_tag(umi_tag):
                            main_missing_umi += 1
                    else:
                        excluded_count += 1

                    # Progress reporting
                    if main_total_processed % report_interval == 0:
                        elapsed = time.time() - start_time
                        rate = main_total_processed / elapsed
                        print(f"  Progress: {main_total_processed:,} reads processed, "
                              f"{main_count:,} written, {excluded_count:,} excluded ({rate:,.0f} reads/sec)",
                              file=sys.stderr, flush=True)

        elapsed = time.time() - start_time
        print(f"Wrote {main_count:,} reads from main BAM in {elapsed:.1f}s", file=sys.stderr)
        if main_missing_cb > 0:
            print(f"  Warning: {main_missing_cb} reads missing {cb_tag} tag", file=sys.stderr)
        if main_missing_umi > 0:
            print(f"  Warning: {main_missing_umi} reads missing {umi_tag} tag", file=sys.stderr)
        if supp_bam:
            print(f"Excluded {excluded_count} reads from main BAM (present in supplemental)", file=sys.stderr)
        print(f"Total reads written: {supp_count + main_count}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description="Extract primary reads from one or two CUDLL BAM files and create CB/UMI table"
    )
    parser.add_argument(
        "--main-bam",
        required=True,
        help="Main CUDLL BAM file"
    )
    parser.add_argument(
        "--supp-bam",
        help="Optional supplemental CUDLL BAM file (takes priority when provided)"
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output tab-delimited file (can be .gz compressed)"
    )
    parser.add_argument(
        "--cb-tag",
        default="CB",
        help="BAM tag for cell barcode (default: CB)"
    )
    parser.add_argument(
        "--umi-tag",
        default="XM",
        help="BAM tag for UMI (default: XM)"
    )

    args = parser.parse_args()

    print(f"Processing CUDLL BAM files...", file=sys.stderr)
    print(f"  Main BAM: {args.main_bam}", file=sys.stderr)
    print(f"  Supp BAM: {args.supp_bam if args.supp_bam else 'none'}", file=sys.stderr)
    print(f"  Output: {args.output}", file=sys.stderr)
    print(f"  Cell Barcode Tag: {args.cb_tag}", file=sys.stderr)
    print(f"  UMI Tag: {args.umi_tag}", file=sys.stderr)

    process_bams(args.main_bam, args.supp_bam, args.output, args.cb_tag, args.umi_tag)

    print("Done!", file=sys.stderr)


if __name__ == "__main__":
    main()
