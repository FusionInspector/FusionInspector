#!/usr/bin/env python3

"""
Extract primary reads from one or two CUDLL BAM files and create a table of read
names with their cell barcodes (CB) and UMIs (XM).
When a supplemental BAM is provided, primary reads from both BAMs are combined
and duplicate read names are collapsed after sorting.
"""

import argparse
import gzip
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pysam


def open_text_writer(path):
    if path.endswith(".gz"):
        return gzip.open(path, "wt")
    return open(path, "wt")


def open_text_reader(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def select_plot_indices(num_points, max_points):
    if num_points <= max_points:
        return list(range(num_points))

    max_log_rank = math.log10(num_points) if num_points > 1 else 1.0
    selected = {0, num_points - 1}

    # Sample ranks across log space so the knee shape remains visible.
    for sample_idx in range(max_points):
        frac = sample_idx / (max_points - 1)
        rank = 10 ** (frac * max_log_rank)
        index = min(num_points - 1, max(0, int(round(rank)) - 1))
        selected.add(index)

    return sorted(selected)


SIZE_UNITS = {
    "": 1,
    "B": 1,
    "K": 1024,
    "KB": 1024,
    "M": 1024 ** 2,
    "MB": 1024 ** 2,
    "G": 1024 ** 3,
    "GB": 1024 ** 3,
    "T": 1024 ** 4,
    "TB": 1024 ** 4,
}


def parse_size_to_bytes(size_text):
    match = re.fullmatch(r"(?i)\s*(\d+)\s*([KMGT]?B?)\s*", size_text)
    if not match:
        raise ValueError(
            f"Invalid size '{size_text}'. Expected forms like 512M, 2G, or 4096."
        )

    value = int(match.group(1))
    unit = match.group(2).upper()
    return value * SIZE_UNITS[unit]


def build_sort_command(sort_exe, sort_threads, sort_memory_per_thread, extra_args):
    total_memory_bytes = sort_threads * parse_size_to_bytes(sort_memory_per_thread)
    return [
        sort_exe,
        "--parallel",
        str(sort_threads),
        "-S",
        str(total_memory_bytes),
        *extra_args,
    ]


class BarcodeUmiAggregator:
    """Spill valid barcode/UMI pairs to disk, sort externally, then summarize."""

    def __init__(
        self,
        summary_file=None,
        knee_plot_file=None,
        max_knee_plot_points=5000,
        sort_exe=None,
        sort_threads=1,
        sort_memory_per_thread="1G",
    ):
        self.summary_file = summary_file
        self.knee_plot_file = knee_plot_file
        self.max_knee_plot_points = max(2, max_knee_plot_points)
        self.enabled = bool(summary_file and knee_plot_file)
        self._sort_exe = sort_exe if self.enabled else None
        self._sort_threads = sort_threads
        self._sort_memory_per_thread = sort_memory_per_thread
        self._pairs_path = None
        self._sorted_path = None
        self._pairs_handle = None

        if self.enabled:
            if self._sort_exe is None:
                raise RuntimeError("GNU sort is required to generate barcode summary outputs")

            pairs_fd, self._pairs_path = tempfile.mkstemp(prefix="barcode_umi_pairs_", suffix=".tsv")
            self._pairs_handle = os.fdopen(pairs_fd, "wt")

    def add(self, cell_barcode, umi):
        if not self.enabled:
            return

        if cell_barcode == "NA" or umi == "NA":
            return

        self._pairs_handle.write(f"{cell_barcode}\t{umi}\n")

    def finalize(self):
        if not self.enabled:
            return

        try:
            self._pairs_handle.close()
            self._sort_pairs()
            barcode_counts = self._write_summary_from_sorted_pairs()
            self._write_knee_plot(barcode_counts)
        finally:
            if self._pairs_handle is not None and not self._pairs_handle.closed:
                self._pairs_handle.close()
            if self._pairs_path and os.path.exists(self._pairs_path):
                os.unlink(self._pairs_path)
            if self._sorted_path and os.path.exists(self._sorted_path):
                os.unlink(self._sorted_path)

    def _sort_pairs(self):
        sorted_fd, self._sorted_path = tempfile.mkstemp(prefix="barcode_umi_pairs_sorted_", suffix=".tsv")
        os.close(sorted_fd)

        subprocess.run(
            build_sort_command(
                self._sort_exe,
                self._sort_threads,
                self._sort_memory_per_thread,
                [
                    "-u",
                    "-t",
                    "\t",
                    "-k1,1",
                    "-k2,2",
                    self._pairs_path,
                    "-o",
                    self._sorted_path,
                ],
            ),
            check=True,
            env={**os.environ, "LC_ALL": "C"},
        )

    def _write_summary_from_sorted_pairs(self):
        counts_fd, counts_path = tempfile.mkstemp(prefix="barcode_umi_counts_", suffix=".tsv")
        os.close(counts_fd)
        sorted_counts_fd, sorted_counts_path = tempfile.mkstemp(
            prefix="barcode_umi_counts_sorted_",
            suffix=".tsv",
        )
        os.close(sorted_counts_fd)

        try:
            current_barcode = None
            current_count = 0

            with open(self._sorted_path, "rt") as input_handle, open(counts_path, "wt") as counts_handle:
                for line in input_handle:
                    cell_barcode, _umi = line.rstrip("\n").split("\t", 1)

                    if cell_barcode != current_barcode:
                        if current_barcode is not None:
                            counts_handle.write(f"{current_barcode}\t{current_count}\n")
                        current_barcode = cell_barcode
                        current_count = 0

                    current_count += 1

                if current_barcode is not None:
                    counts_handle.write(f"{current_barcode}\t{current_count}\n")

            subprocess.run(
                build_sort_command(
                    self._sort_exe,
                    self._sort_threads,
                    self._sort_memory_per_thread,
                    [
                        "-t",
                        "\t",
                        "-k2,2nr",
                        "-k1,1",
                        counts_path,
                        "-o",
                        sorted_counts_path,
                    ],
                ),
                check=True,
                env={**os.environ, "LC_ALL": "C"},
            )

            barcode_counts = []
            with open(sorted_counts_path, "rt") as summary_fh:
                for line in summary_fh:
                    cell_barcode, umi_count = line.rstrip("\n").split("\t", 1)
                    barcode_counts.append((cell_barcode, int(umi_count)))

            with open_text_writer(self.summary_file) as summary_fh:
                summary_fh.write("cell_barcode\tumi_count\n")
                for cell_barcode, umi_count in barcode_counts:
                    summary_fh.write(f"{cell_barcode}\t{umi_count}\n")

            return barcode_counts
        finally:
            if os.path.exists(counts_path):
                os.unlink(counts_path)
            if os.path.exists(sorted_counts_path):
                os.unlink(sorted_counts_path)

    def _write_knee_plot(self, barcode_counts):
        fig, ax = plt.subplots(figsize=(7, 5))

        if barcode_counts:
            umi_counts = [count for _, count in barcode_counts]
            plot_indices = select_plot_indices(len(umi_counts), self.max_knee_plot_points)
            ranks = [index + 1 for index in plot_indices]
            umi_counts = [umi_counts[index] for index in plot_indices]
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


def get_cb_umi_fields(read, cb_tag, umi_tag):
    try:
        cell_barcode = read.get_tag(cb_tag)
    except KeyError:
        cell_barcode = "NA"

    try:
        umi = read.get_tag(umi_tag)
    except KeyError:
        umi = "NA"

    return read.query_name, cell_barcode, umi


def append_primary_reads_to_combined_file(combined_handle, bam_path, cb_tag, umi_tag, label):
    primary_count = 0
    report_interval = 1000000  # Report every 1M reads
    start_time = time.time()

    print(f"Reading primary reads from {label} BAM: {bam_path}", file=sys.stderr)

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam:
            if read.is_secondary or read.is_supplementary:
                continue

            read_name, cell_barcode, umi = get_cb_umi_fields(read, cb_tag, umi_tag)
            combined_handle.write(f"{read_name}\t{cell_barcode}\t{umi}\n")
            primary_count += 1

            if primary_count % report_interval == 0:
                elapsed = time.time() - start_time
                rate = primary_count / elapsed
                print(
                    f"  Progress: {primary_count:,} primary reads written from {label} BAM ({rate:,.0f} reads/sec)",
                    file=sys.stderr,
                    flush=True,
                )

    elapsed = time.time() - start_time
    print(
        f"Collected {primary_count:,} primary reads from {label} BAM in {elapsed:.1f}s",
        file=sys.stderr,
    )
    return primary_count


def sort_combined_records(sort_exe, sort_threads, sort_memory_per_thread, combined_path):
    sorted_fd, sorted_path = tempfile.mkstemp(prefix="cb_umi_records_sorted_", suffix=".tsv")
    os.close(sorted_fd)

    subprocess.run(
        build_sort_command(
            sort_exe,
            sort_threads,
            sort_memory_per_thread,
            [
                "-t",
                "\t",
                "-k1,1",
                "-k2,2",
                "-k3,3",
                combined_path,
                "-o",
                sorted_path,
            ],
        ),
        check=True,
        env={**os.environ, "LC_ALL": "C"},
    )

    return sorted_path


def process_bams(
    main_bam,
    supp_bam,
    output_file,
    cb_tag,
    umi_tag,
    summary_file=None,
    knee_plot_file=None,
    max_knee_plot_points=5000,
    sort_threads=1,
    sort_memory_per_thread="1G",
):
    """
    Process one or two BAM files and create a CB/UMI table.
    When supp_bam is provided, records are merged and deduplicated by read name.
    """
    sort_exe = shutil.which("sort")
    if sort_exe is None:
        raise RuntimeError("GNU sort is required")

    if sort_threads < 1:
        raise ValueError("--sort-threads must be at least 1")
    parse_size_to_bytes(sort_memory_per_thread)

    aggregator = BarcodeUmiAggregator(
        summary_file,
        knee_plot_file,
        max_knee_plot_points,
        sort_exe=sort_exe,
        sort_threads=sort_threads,
        sort_memory_per_thread=sort_memory_per_thread,
    )

    combined_fd, combined_path = tempfile.mkstemp(prefix="cb_umi_records_", suffix=".tsv")
    sorted_path = None

    try:
        with os.fdopen(combined_fd, "wt") as combined_handle:
            if supp_bam:
                append_primary_reads_to_combined_file(
                    combined_handle,
                    supp_bam,
                    cb_tag,
                    umi_tag,
                    "supplemental",
                )
            else:
                print("No supplemental BAM provided; processing only the main BAM", file=sys.stderr)

            append_primary_reads_to_combined_file(
                combined_handle,
                main_bam,
                cb_tag,
                umi_tag,
                "main",
            )

        print("Sorting combined read-level records by read name", file=sys.stderr)
        sorted_path = sort_combined_records(
            sort_exe,
            sort_threads,
            sort_memory_per_thread,
            combined_path,
        )

        opener = gzip.open if output_file.endswith(".gz") else open
        retained_count = 0
        duplicate_count = 0
        missing_cb_count = 0
        missing_umi_count = 0
        previous_read_name = None
        start_time = time.time()
        report_interval = 1000000  # Report every 1M retained reads

        with open(sorted_path, "rt") as sorted_handle, opener(output_file, "wt") as outfile:
            outfile.write("read_name\tcell_barcode\tUMI\n")

            for line in sorted_handle:
                read_name, cell_barcode, umi = line.rstrip("\n").split("\t")
                if read_name == previous_read_name:
                    duplicate_count += 1
                    continue

                previous_read_name = read_name
                outfile.write(f"{read_name}\t{cell_barcode}\t{umi}\n")
                aggregator.add(cell_barcode, umi)
                retained_count += 1

                if cell_barcode == "NA":
                    missing_cb_count += 1
                if umi == "NA":
                    missing_umi_count += 1

                if retained_count % report_interval == 0:
                    elapsed = time.time() - start_time
                    rate = retained_count / elapsed
                    print(
                        f"  Progress: {retained_count:,} unique read names written ({rate:,.0f} reads/sec)",
                        file=sys.stderr,
                        flush=True,
                    )

        elapsed = time.time() - start_time
        print(f"Wrote {retained_count:,} unique read-level records in {elapsed:.1f}s", file=sys.stderr)
        if duplicate_count > 0:
            print(f"Discarded {duplicate_count:,} duplicate records by read name", file=sys.stderr)
        if missing_cb_count > 0:
            print(f"  Warning: {missing_cb_count} retained reads missing {cb_tag} tag", file=sys.stderr)
        if missing_umi_count > 0:
            print(f"  Warning: {missing_umi_count} retained reads missing {umi_tag} tag", file=sys.stderr)
    finally:
        if os.path.exists(combined_path):
            os.unlink(combined_path)
        if sorted_path and os.path.exists(sorted_path):
            os.unlink(sorted_path)

    if aggregator.enabled:
        print(f"Sorting barcode/UMI pairs for summary generation", file=sys.stderr)
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
        help="Optional supplemental CUDLL BAM file to merge with the main BAM before deduplicating by read name"
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output tab-delimited file (can be .gz compressed)"
    )
    parser.add_argument(
        "--barcode-umi-counts",
        help="Output TSV containing unique cell barcodes and their unique UMI counts"
    )
    parser.add_argument(
        "--knee-plot-pdf",
        help="Output PDF for the barcode-rank knee plot"
    )
    parser.add_argument(
        "--max-knee-plot-points",
        type=int,
        default=5000,
        help="Maximum number of barcode ranks to render in the knee plot (default: 5000)"
    )
    parser.add_argument(
        "--sort-threads",
        type=int,
        default=1,
        help="Number of threads to use for each external sort invocation (default: 1)"
    )
    parser.add_argument(
        "--sort-memory-per-thread",
        default="1G",
        help="Approximate RAM budget per sort thread, for example 512M or 2G (default: 1G)"
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

    if bool(args.barcode_umi_counts) != bool(args.knee_plot_pdf):
        parser.error("--barcode-umi-counts and --knee-plot-pdf must be provided together")

    print("Processing CUDLL BAM files...", file=sys.stderr)
    print(f"  Main BAM: {args.main_bam}", file=sys.stderr)
    print(f"  Supp BAM: {args.supp_bam if args.supp_bam else 'none'}", file=sys.stderr)
    print(f"  Output: {args.output}", file=sys.stderr)
    print(f"  Barcode UMI Counts: {args.barcode_umi_counts if args.barcode_umi_counts else 'disabled'}", file=sys.stderr)
    print(f"  Knee Plot PDF: {args.knee_plot_pdf if args.knee_plot_pdf else 'disabled'}", file=sys.stderr)
    print(f"  Cell Barcode Tag: {args.cb_tag}", file=sys.stderr)
    print(f"  UMI Tag: {args.umi_tag}", file=sys.stderr)
    print(f"  Max Knee Plot Points: {args.max_knee_plot_points}", file=sys.stderr)
    print(f"  Sort Threads: {args.sort_threads}", file=sys.stderr)
    print(f"  Sort Memory Per Thread: {args.sort_memory_per_thread}", file=sys.stderr)

    process_bams(
        args.main_bam,
        args.supp_bam,
        args.output,
        args.cb_tag,
        args.umi_tag,
        args.barcode_umi_counts,
        args.knee_plot_pdf,
        args.max_knee_plot_points,
        args.sort_threads,
        args.sort_memory_per_thread,
    )

    print("Done!", file=sys.stderr)


if __name__ == "__main__":
    main()
