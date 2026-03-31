#!/usr/bin/env python3

import csv
import gzip
import os
import sys


EXPECTED_CB_UMI_ROWS = [
    ("readA", "CELL_A", "UMI_1"),
    ("readB", "CELL_A", "UMI_1"),
    ("readC", "CELL_A", "UMI_2"),
    ("readD", "CELL_B", "UMI_3"),
    ("readE", "CELL_C", "NA"),
    ("readF", "NA", "UMI_4"),
    ("shared1", "CELL_A", "UMI_7"),
    ("shared2", "CELL_B", "UMI_3"),
    ("supp_only1", "CELL_B", "UMI_8"),
    ("supp_only2", "CELL_E", "UMI_10"),
]

EXPECTED_BARCODE_COUNTS = [
    ("CELL_A", "3"),
    ("CELL_B", "2"),
    ("CELL_E", "1"),
]


def fail(message):
    print(f"ERROR: {message}", file=sys.stderr)
    sys.exit(1)


def read_cb_umi_table(path):
    with gzip.open(path, "rt") as handle:
        reader = csv.reader(handle, delimiter="\t")
        rows = list(reader)

    if not rows:
        fail("CB/UMI output is empty")

    header = rows[0]
    if header != ["read_name", "cell_barcode", "UMI"]:
        fail(f"Unexpected CB/UMI header: {header}")

    return [tuple(row) for row in rows[1:]]


def read_barcode_counts(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as handle:
        reader = csv.reader(handle, delimiter="\t")
        rows = list(reader)

    if not rows:
        fail("barcode count output is empty")

    header = rows[0]
    if header != ["cell_barcode", "umi_count"]:
        fail(f"Unexpected barcode-count header: {header}")

    return [tuple(row) for row in rows[1:]]


def main():
    if len(sys.argv) not in (2, 4):
        sys.exit(
            f"usage: {sys.argv[0]} <cb_umi.tsv.gz> [<barcode_counts.tsv> <knee.pdf>]"
        )

    cb_umi_path = sys.argv[1]

    cb_umi_rows = read_cb_umi_table(cb_umi_path)
    if cb_umi_rows != EXPECTED_CB_UMI_ROWS:
        fail(f"CB/UMI rows differ.\nExpected: {EXPECTED_CB_UMI_ROWS}\nObserved: {cb_umi_rows}")

    if len(sys.argv) == 2:
        print("Primary CB/UMI output validation passed.")
        return

    barcode_counts_path, knee_pdf_path = sys.argv[2:]

    barcode_counts = read_barcode_counts(barcode_counts_path)
    if barcode_counts != EXPECTED_BARCODE_COUNTS:
        fail(
            "barcode-count rows differ.\n"
            f"Expected: {EXPECTED_BARCODE_COUNTS}\nObserved: {barcode_counts}"
        )

    if not os.path.exists(knee_pdf_path):
        fail(f"knee plot PDF not found: {knee_pdf_path}")
    if os.path.getsize(knee_pdf_path) == 0:
        fail(f"knee plot PDF is empty: {knee_pdf_path}")

    print("Execution test passed.")


if __name__ == "__main__":
    main()
