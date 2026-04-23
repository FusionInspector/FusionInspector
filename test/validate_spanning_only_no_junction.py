#!/usr/bin/env python3

import csv
import os
import subprocess
import sys


def main():
    if len(sys.argv) != 3:
        sys.exit(
            "usage: validate_spanning_only_no_junction.py output_dir expected_fusion"
        )

    output_dir = sys.argv[1]
    expected_fusion = sys.argv[2]
    fusion_tsv = os.path.join(output_dir, "finspector.FusionInspector.fusions.tsv")
    junction_bam = os.path.join(output_dir, "IGV_inputs", "finspector.junction_reads.bam")
    spanning_bam = os.path.join(output_dir, "IGV_inputs", "finspector.spanning_reads.bam")

    row = find_fusion_row(fusion_tsv, expected_fusion)
    junction_count = int(row["JunctionReadCount"])
    spanning_count = int(row["SpanningFragCount"])

    if junction_count != 0:
        sys.exit(f"Error, expected zero junction reads but found {junction_count}")
    if spanning_count <= 0:
        sys.exit(f"Error, expected spanning support but found {spanning_count}")

    samtools_quickcheck(junction_bam)
    samtools_quickcheck(spanning_bam)

    junction_bam_records = samtools_view_count(junction_bam)
    spanning_bam_records = samtools_view_count(spanning_bam)

    if junction_bam_records != 0:
        sys.exit(
            f"Error, expected empty junction BAM but found {junction_bam_records} records"
        )
    if spanning_bam_records <= 0:
        sys.exit(
            f"Error, expected non-empty spanning BAM but found {spanning_bam_records} records"
        )

    print(
        "Validated spanning-only test: "
        f"{expected_fusion} has J={junction_count}, S={spanning_count}, "
        f"junction_bam_records={junction_bam_records}, "
        f"spanning_bam_records={spanning_bam_records}"
    )


def find_fusion_row(fusion_tsv, expected_fusion):
    with open(fusion_tsv, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        matches = [row for row in reader if row["#FusionName"] == expected_fusion]

    if not matches:
        sys.exit(f"Error, fusion {expected_fusion} not found in {fusion_tsv}")

    matches.sort(
        key=lambda row: (int(row["JunctionReadCount"]), -int(row["SpanningFragCount"]))
    )
    return matches[0]


def samtools_quickcheck(bam_file):
    subprocess.run(["samtools", "quickcheck", bam_file], check=True)


def samtools_view_count(bam_file):
    result = subprocess.run(
        ["samtools", "view", "-c", bam_file],
        check=True,
        text=True,
        stdout=subprocess.PIPE,
    )
    return int(result.stdout.strip())


if __name__ == "__main__":
    main()
