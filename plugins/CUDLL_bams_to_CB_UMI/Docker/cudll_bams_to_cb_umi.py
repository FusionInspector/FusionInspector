#!/usr/bin/env python3

"""
Extract primary reads from one or two CUDLL BAM files and create a table of read
names with their cell barcodes (CB) and UMIs (XM).
When provided, reads from the supplemental BAM take priority; reads with the
same name in the main BAM are excluded.
"""

import argparse
import gzip
import os
import sqlite3
import sys
import tempfile
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pysam


class BarcodeUmiAggregator:
    """Track unique barcode/UMI pairs on disk and emit summary outputs."""

    def __init__(self, summary_file, knee_plot_file):
        self.summary_file = summary_file
        self.knee_plot_file = knee_plot_file
        self.db_fd, self.db_path = tempfile.mkstemp(prefix="barcode_umi_", suffix=".sqlite3")
        os.close(self.db_fd)
        self.conn = sqlite3.connect(self.db_path)
        self.conn.execute("PRAGMA journal_mode = OFF")
        self.conn.execute("PRAGMA synchronous = OFF")
        self.conn.execute("PRAGMA temp_store = MEMORY")
        self.conn.execute(
            """
            CREATE TABLE barcode_umi (
                cell_barcode TEXT NOT NULL,
                umi TEXT NOT NULL,
                PRIMARY KEY (cell_barcode, umi)
            )
            """
        )
        self._pending_records = []
        self._batch_size = 10000

    def add(self, cell_barcode, umi):
        if cell_barcode == "NA" or umi == "NA":
            return

        self._pending_records.append((cell_barcode, umi))
        if len(self._pending_records) >= self._batch_size:
            self.flush()

    def flush(self):
        if not self._pending_records:
            return

        self.conn.executemany(
            "INSERT OR IGNORE INTO barcode_umi(cell_barcode, umi) VALUES (?, ?)",
            self._pending_records,
        )
        self.conn.commit()
        self._pending_records.clear()

    def finalize(self):
        self.flush()

        cursor = self.conn.execute(
            """
            SELECT cell_barcode, COUNT(*) AS umi_count
            FROM barcode_umi
            GROUP BY cell_barcode
            ORDER BY umi_count DESC, cell_barcode
            """
        )
        barcode_counts = cursor.fetchall()

        with open(self.summary_file, "wt") as summary_fh:
            summary_fh.write("cell_barcode\tumi_count\n")
            for cell_barcode, umi_count in barcode_counts:
                summary_fh.write(f"{cell_barcode}\t{umi_count}\n")

        self._write_knee_plot(barcode_counts)

        self.conn.close()
        os.unlink(self.db_path)

    def _write_knee_plot(self, barcode_counts):
        fig, ax = plt.subplots(figsize=(7, 5))

        if barcode_counts:
            ranks = list(range(1, len(barcode_counts) + 1))
            umi_counts = [count for _, count in barcode_counts]
            ax.loglog(ranks, umi_counts, linewidth=1.5)
            ax.set_xlim(left=1)
            ax.set_ylim(bottom=1)
        else:
            ax.text(
                0.5,
                0.5,
                "No valid barcode/UMI pairs found",
                ha="center",
                va="center",
                transform=ax.transAxes,
            )

        ax.set_title("Barcode UMI Knee Plot")
        ax.set_xlabel("Barcode rank")
        ax.set_ylabel("Unique UMI count")
        ax.grid(True, which="both", linestyle=":", linewidth=0.5)

        fig.tight_layout()
        fig.savefig(self.knee_plot_file)
        plt.close(fig)


def extract_read_names(bam_path):
    """Extract all primary read names from a BAM file."""
    read_names = set()
    start_time = time.time()
    total_reads = 0
    report_interval = 1000000  # Report every 1M reads

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam:
            total_reads += 1
            if not read.is_secondary and not read.is_supplementary:
                read_names.add(read.query_name)

            if total_reads % report_interval == 0:
                elapsed = time.time() - start_time
                rate = total_reads / elapsed
                print(
                    f"  Progress: {total_reads:,} reads processed ({rate:,.0f} reads/sec)",
                    file=sys.stderr,
                    flush=True,
                )

    elapsed = time.time() - start_time
    print(
        f"  Completed: {total_reads:,} total reads in {elapsed:.1f}s",
        file=sys.stderr,
        flush=True,
    )
    return read_names


def write_cb_umi_record(outfile, read, cb_tag, umi_tag, aggregator):
    """Write a single read's CB/UMI info to the output file."""
    read_name = read.query_name

    try:
        cell_barcode = read.get_tag(cb_tag)
    except KeyError:
        cell_barcode = "NA"

    try:
        umi = read.get_tag(umi_tag)
    except KeyError:
        umi = "NA"

    outfile.write(f"{read_name}\t{cell_barcode}\t{umi}\n")
    aggregator.add(cell_barcode, umi)


def process_bams(main_bam, supp_bam, output_file, cb_tag, umi_tag, summary_file, knee_plot_file):
    """
    Process one or two BAM files and create a CB/UMI table.
    Priority is given to reads from supp_bam when provided.
    """
    supp_read_names = set()
    supp_count = 0
    report_interval = 1000000  # Report every 1M reads
    aggregator = BarcodeUmiAggregator(summary_file, knee_plot_file)

    if supp_bam:
        print(f"Extracting read names from supplemental BAM: {supp_bam}", file=sys.stderr)
        supp_read_names = extract_read_names(supp_bam)
        print(f"Found {len(supp_read_names)} unique read names in supplemental BAM", file=sys.stderr)
    else:
        print("No supplemental BAM provided; processing only the main BAM", file=sys.stderr)

    opener = gzip.open if output_file.endswith(".gz") else open
    with opener(output_file, "wt") as outfile:
        outfile.write("read_name\tcell_barcode\tUMI\n")

        if supp_bam:
            print("Writing reads from supplemental BAM...", file=sys.stderr)
            supp_missing_cb = 0
            supp_missing_umi = 0
            start_time = time.time()

            with pysam.AlignmentFile(supp_bam, "rb") as bam:
                for read in bam:
                    if not read.is_secondary and not read.is_supplementary:
                        write_cb_umi_record(outfile, read, cb_tag, umi_tag, aggregator)
                        supp_count += 1
                        if not read.has_tag(cb_tag):
                            supp_missing_cb += 1
                        if not read.has_tag(umi_tag):
                            supp_missing_umi += 1

                        if supp_count % report_interval == 0:
                            elapsed = time.time() - start_time
                            rate = supp_count / elapsed
                            print(
                                f"  Progress: {supp_count:,} reads written ({rate:,.0f} reads/sec)",
                                file=sys.stderr,
                                flush=True,
                            )

            elapsed = time.time() - start_time
            print(f"Wrote {supp_count:,} reads from supplemental BAM in {elapsed:.1f}s", file=sys.stderr)
            if supp_missing_cb > 0:
                print(f"  Warning: {supp_missing_cb} reads missing {cb_tag} tag", file=sys.stderr)
            if supp_missing_umi > 0:
                print(f"  Warning: {supp_missing_umi} reads missing {umi_tag} tag", file=sys.stderr)

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
                        write_cb_umi_record(outfile, read, cb_tag, umi_tag, aggregator)
                        main_count += 1
                        if not read.has_tag(cb_tag):
                            main_missing_cb += 1
                        if not read.has_tag(umi_tag):
                            main_missing_umi += 1
                    else:
                        excluded_count += 1

                    if main_total_processed % report_interval == 0:
                        elapsed = time.time() - start_time
                        rate = main_total_processed / elapsed
                        print(
                            f"  Progress: {main_total_processed:,} reads processed, "
                            f"{main_count:,} written, {excluded_count:,} excluded ({rate:,.0f} reads/sec)",
                            file=sys.stderr,
                            flush=True,
                        )

        elapsed = time.time() - start_time
        print(f"Wrote {main_count:,} reads from main BAM in {elapsed:.1f}s", file=sys.stderr)
        if main_missing_cb > 0:
            print(f"  Warning: {main_missing_cb} reads missing {cb_tag} tag", file=sys.stderr)
        if main_missing_umi > 0:
            print(f"  Warning: {main_missing_umi} reads missing {umi_tag} tag", file=sys.stderr)
        if supp_bam:
            print(f"Excluded {excluded_count} reads from main BAM (present in supplemental)", file=sys.stderr)
        print(f"Total reads written: {supp_count + main_count}", file=sys.stderr)

    print(f"Writing barcode UMI summary: {summary_file}", file=sys.stderr)
    print(f"Writing barcode knee plot PDF: {knee_plot_file}", file=sys.stderr)
    aggregator.finalize()


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
        "--barcode-umi-counts",
        required=True,
        help="Output TSV containing unique cell barcodes and their unique UMI counts"
    )
    parser.add_argument(
        "--knee-plot-pdf",
        required=True,
        help="Output PDF for the barcode-rank knee plot"
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

    print("Processing CUDLL BAM files...", file=sys.stderr)
    print(f"  Main BAM: {args.main_bam}", file=sys.stderr)
    print(f"  Supp BAM: {args.supp_bam if args.supp_bam else 'none'}", file=sys.stderr)
    print(f"  Output: {args.output}", file=sys.stderr)
    print(f"  Barcode UMI Counts: {args.barcode_umi_counts}", file=sys.stderr)
    print(f"  Knee Plot PDF: {args.knee_plot_pdf}", file=sys.stderr)
    print(f"  Cell Barcode Tag: {args.cb_tag}", file=sys.stderr)
    print(f"  UMI Tag: {args.umi_tag}", file=sys.stderr)

    process_bams(
        args.main_bam,
        args.supp_bam,
        args.output,
        args.cb_tag,
        args.umi_tag,
        args.barcode_umi_counts,
        args.knee_plot_pdf,
    )

    print("Done!", file=sys.stderr)


if __name__ == "__main__":
    main()
