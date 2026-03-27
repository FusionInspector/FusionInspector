# Execution Test

This directory contains a deterministic local execution test for the `CUDLL_bams_to_CB_UMI` plugin.

## What it does

`run_execution_test.sh` performs three steps:

1. Generates small synthetic main and supplemental BAMs with:
   - overlapping read names across BAMs
   - duplicate UMIs within the same barcode
   - missing CB or UMI tags
   - secondary and supplementary alignments that should be ignored
2. Runs `Docker/cudll_bams_to_cb_umi.py` on those BAMs.
3. Validates:
   - exact read-level CB/UMI table content
   - exact barcode-level unique UMI counts
   - existence of a non-empty knee-plot PDF

## Run

```bash
cd plugins/CUDLL_bams_to_CB_UMI
./testing/run_execution_test.sh
```

## Layout

- `bin/generate_simulated_bams.py`: creates the synthetic BAM inputs
- `bin/validate_execution_test.py`: checks the plugin outputs against expected values
- `tmp/`: regenerated synthetic BAMs
- `output/`: regenerated plugin outputs
