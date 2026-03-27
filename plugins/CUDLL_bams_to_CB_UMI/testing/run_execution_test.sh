#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PLUGIN_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
TMP_DIR="${PLUGIN_DIR}/testing/tmp"
OUT_DIR="${PLUGIN_DIR}/testing/output"

rm -rf "${TMP_DIR}" "${OUT_DIR}"
mkdir -p "${TMP_DIR}" "${OUT_DIR}"

python3 "${PLUGIN_DIR}/testing/bin/generate_simulated_bams.py" "${TMP_DIR}" >/dev/null

python3 "${PLUGIN_DIR}/Docker/cudll_bams_to_cb_umi.py" \
    --main-bam "${TMP_DIR}/sim.main.bam" \
    --supp-bam "${TMP_DIR}/sim.supp.bam" \
    --output "${OUT_DIR}/sim.cb_umi.tsv.gz" \
    --barcode-umi-counts "${OUT_DIR}/sim.barcode_umi_counts.tsv" \
    --knee-plot-pdf "${OUT_DIR}/sim.barcode_knee_plot.pdf" \
    --cb-tag CB \
    --umi-tag XM

python3 "${PLUGIN_DIR}/testing/bin/validate_execution_test.py" \
    "${OUT_DIR}/sim.cb_umi.tsv.gz" \
    "${OUT_DIR}/sim.barcode_umi_counts.tsv" \
    "${OUT_DIR}/sim.barcode_knee_plot.pdf"
