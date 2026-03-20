#!/usr/bin/env python3

import argparse
import csv
import sys


def main():
    csv.field_size_limit(sys.maxsize)

    parser = argparse.ArgumentParser(
        description="Filter long-read fusions using ctat-LR-fusion-like abundance criteria."
    )
    parser.add_argument("--fusion_preds", required=True, help="Fusion prediction summary TSV")
    parser.add_argument(
        "--min_LR_reads",
        type=int,
        default=1,
        help="Minimum long-read support for reference-spliced events",
    )
    parser.add_argument(
        "--min_LR_novel_reads",
        type=int,
        default=2,
        help="Minimum long-read support for non-reference-spliced events",
    )
    args = parser.parse_args()

    with open(args.fusion_preds) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None:
            return
        writer = csv.DictWriter(
            sys.stdout, fieldnames=reader.fieldnames, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()

        for row in reader:
            j = int(float(row.get("JunctionReadCount", "0") or 0))
            splice_type = row.get("SpliceType", "")
            if splice_type == "ONLY_REF_SPLICE":
                keep = j >= args.min_LR_reads
            else:
                keep = j >= args.min_LR_novel_reads

            if keep:
                writer.writerow(row)


if __name__ == "__main__":
    main()
