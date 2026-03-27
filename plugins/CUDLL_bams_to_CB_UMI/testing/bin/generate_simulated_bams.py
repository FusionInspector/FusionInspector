#!/usr/bin/env python3

import os
import sys

import pysam


HEADER = {
    "HD": {"VN": "1.6", "SO": "unsorted"},
    "SQ": [{"SN": "chrTest", "LN": 100000}],
}


def make_read(name, start, cb=None, umi=None, *, secondary=False, supplementary=False):
    read = pysam.AlignedSegment()
    read.query_name = name
    read.query_sequence = "A" * 50
    read.flag = 0
    if secondary:
        read.flag |= 0x100
    if supplementary:
        read.flag |= 0x800
    read.reference_id = 0
    read.reference_start = start
    read.mapping_quality = 60
    read.cigar = ((0, 50),)
    read.query_qualities = pysam.qualitystring_to_array("I" * 50)
    if cb is not None:
        read.set_tag("CB", cb)
    if umi is not None:
        read.set_tag("XM", umi)
    return read


def write_bam(path, reads):
    with pysam.AlignmentFile(path, "wb", header=HEADER) as bam:
        for read in reads:
            bam.write(read)

    pysam.index(path)


def main():
    if len(sys.argv) != 2:
        sys.exit(f"usage: {sys.argv[0]} <output_dir>")

    outdir = sys.argv[1]
    os.makedirs(outdir, exist_ok=True)

    main_bam = os.path.join(outdir, "sim.main.bam")
    supp_bam = os.path.join(outdir, "sim.supp.bam")

    main_reads = [
        make_read("readA", 100, "CELL_A", "UMI_1"),
        make_read("readB", 200, "CELL_A", "UMI_1"),
        make_read("readC", 300, "CELL_A", "UMI_2"),
        make_read("readD", 400, "CELL_B", "UMI_3"),
        make_read("readE", 500, "CELL_C"),
        make_read("readF", 600, umi="UMI_4"),
        make_read("readG", 700, "CELL_D", "UMI_5", supplementary=True),
        make_read("readH", 800, "CELL_D", "UMI_6", secondary=True),
        make_read("shared1", 900, "CELL_X", "UMI_9"),
        make_read("shared2", 1000, "CELL_Z", "UMI_11"),
    ]

    supp_reads = [
        make_read("shared1", 1100, "CELL_A", "UMI_7"),
        make_read("shared2", 1200, "CELL_B", "UMI_3"),
        make_read("supp_only1", 1300, "CELL_B", "UMI_8"),
        make_read("supp_only2", 1400, "CELL_E", "UMI_10"),
        make_read("supp_secondary", 1500, "CELL_Q", "UMI_12", secondary=True),
    ]

    write_bam(main_bam, main_reads)
    write_bam(supp_bam, supp_reads)

    print(main_bam)
    print(supp_bam)


if __name__ == "__main__":
    main()
