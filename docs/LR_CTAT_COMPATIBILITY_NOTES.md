# FusionInspector Long-Read ctat-LR-fusion Compatibility Notes

## Purpose

This document records the long-read-specific FusionInspector changes made to better match the evidence-capture behavior of `ctat-LR-fusion` while preserving FusionInspector downstream characterization.

The immediate objective was:

- if a fusion is present in the `ctat-LR-fusion` long-read candidate list supplied to FusionInspector,
- and the same fusion contigs / remapped long-read alignments are being evaluated,
- then FusionInspector should retain that fusion instead of dropping it during the original short-read-oriented evidence filtering stages.


## Background

The original FusionInspector long-read path still used the standard FI junction/spanning extraction logic:

- `util/get_fusion_JUNCTION_reads_from_fusion_contig_bam.pl`
- `util/get_fusion_SPANNING_reads_from_bam.from_chim_summary.pl`
- `util/filter_fusions_by_frag_thresholds.pl`

That path is tuned for short-read evidence semantics:

- junction reads
- spanning fragments
- anchor-based filtering
- LDAS gating
- EM-adjusted `J` / `S`

`ctat-LR-fusion` phase 2 instead uses a different long-read evidence model:

- remap candidate long reads to FI-like contigs
- convert remapped long-read alignments to GFF3, allowing non-primary alignments
- count long reads crossing both fusion partners with sufficient exon-overlap support
- filter mostly on long-read count and splice type

This mismatch was the main reason many long-read-supported fusions were being dropped by FusionInspector.


## Files Modified

### Main driver

- `FusionInspector`

### New long-read utilities

- `util/LR_SAM_to_gff3.pl`
- `util/LR_capture_fusion_support_from_gff3.pl`
- `util/LR_filter_fusions_by_evidence_abundance.py`

### Supporting bug fix

- `FusionFilter/util/promiscuity_filter.pl`


## Current Long-Read Flow

When `--read_type long` is used, FusionInspector now branches into a long-read-specific support-capture path.

### 1. Alignment conversion

`util/LR_SAM_to_gff3.pl`

- reads the remapped long-read BAM
- allows non-primary alignments
- converts contig alignments into GFF3-style `match` features
- computes percent identity from `NM` and CIGAR

This is modeled on the analogous conversion used by `ctat-LR-fusion`.

### 2. Long-read support capture

`util/LR_capture_fusion_support_from_gff3.pl`

- parses the FI contig GTF
- tracks original genomic coordinates and splice-site mappings
- parses long-read GFF3 alignments
- excludes seq-similar regions
- requires:
  - support on both fusion partners
  - exon-overlapping alignment support on both sides
  - minimum transcript overlap per side
- snaps local breakpoints back to nearby annotated splice-compatible positions

Output is FI-style summary rows with columns including:

- `JunctionReadCount`
- `SpanningFragCount`
- `est_J`
- `est_S`
- `LargeAnchorSupport`
- `JunctionReads`
- `NumCounterFusionLeft`
- `NumCounterFusionRight`
- `FAR_left`
- `FAR_right`

### 3. Long-read abundance filtering

`util/LR_filter_fusions_by_evidence_abundance.py`

This replaces the short-read fragment-threshold filter for `--read_type long`.

Current thresholds:

- `ONLY_REF_SPLICE`: keep if `JunctionReadCount >= min_LR_reads`
- non-reference splice: keep if `JunctionReadCount >= min_LR_novel_reads`

Current defaults:

- `--min_LR_reads 1`
- `--min_LR_novel_reads 2`

### 4. EM bypass for long reads

For `--read_type long`, FusionInspector now skips both EM passes.

Reason:

- every long read is treated as breakpoint-specific evidence here
- there is no need to re-estimate `J` and `S`

Instead:

- `est_J = JunctionReadCount`
- `est_S = 0.00`

These values are written directly by `util/LR_capture_fusion_support_from_gff3.pl`.

### 5. Downstream characterization retained

After the LR support table is produced, the existing downstream FI steps still run:

- splice annotation
- blast filter
- promiscuity filter
- FFPM incorporation
- microhomology annotation
- FusionAnnotator annotation
- final report generation


## Output Conventions For Long Reads

To minimize downstream code changes and preserve FI output shape:

- `JunctionReadCount` = number of long reads
- `SpanningFragCount` = `0`
- `est_J` = `JunctionReadCount`
- `est_S` = `0.00`
- `LargeAnchorSupport` = `YES`
- `JunctionReads` = comma-delimited long-read accessions
- `SpanningFrags` = `.`


## Counter-Evidence Status

### Implemented

`util/LR_capture_fusion_support_from_gff3.pl` now includes a first-pass long-read counter-evidence calculation.

Current rule:

- left contrary support:
  - read spans the left breakpoint
  - read does not extend into the right gene portion
- right contrary support:
  - read spans the right breakpoint
  - read starts after the left gene portion

This gives non-zero counter support for some long-read events and avoids all-FAR values being purely pseudocount-derived.

### Limitations

This is still a simplified approximation of true long-read counter-support logic.

It does **not yet** perform a full transcript-aware “normal isoform compatibility” analysis analogous to what might ideally be done for long reads.

As a result:

- many events still have `NumCounterFusionLeft = 0` and `NumCounterFusionRight = 0`
- FAR values are still often dominated by the pseudocount


## Validation Performed

Primary validation dataset:

- `/home/unix/bhaas/projects/Melanoma/FUSION-reads-only`

Key comparison inputs:

- `ctat_LR_fusion_outdir/ctat-LR-fusion.fusion_predictions.tsv`
- `fusions.list`
- `FI-min70pid`

Validation outcomes:

### Original FI long-read behavior

Directory:

- `/home/unix/bhaas/projects/Melanoma/FUSION-reads-only/FI-min70pid`

Observed final fusion names:

- 10 unique names

### LR-compatible behavior

Directory:

- `/home/unix/bhaas/projects/Melanoma/FUSION-reads-only/FI-min70pid-lrcompat-noEM-counters`

Observed final fusion names:

- 27 unique names

This matched all fusion names supplied to FusionInspector in:

- `/home/unix/bhaas/projects/Melanoma/FUSION-reads-only/fusions.list`

Important note:

- `ctat-LR-fusion` had one additional fusion, `FRS2--TMEM178B`
- it was **not** present in `fusions.list`
- therefore FusionInspector never built a contig for it
- its absence from final FI output is expected for this test


## Bug Fixes Exposed During LR Adaptation

### Promiscuity filter zero-value bug

File:

- `FusionFilter/util/promiscuity_filter.pl`

Issue:

- code used `or confess` on `sum_JS`
- rows with `sum_JS == 0` were treated as failure
- this caused the long-read-compatible run to die during promiscuity filtering

Fix:

- changed the check to `unless defined $sum_support`

This is a general robustness fix, not just a long-read fix.


## Known Remaining Gaps

### 1. Counter-evidence needs a stronger LR model

The current LR counter-read logic is intentionally simple and should be revisited.

Recommended future direction:

- evaluate long-read counter-support using transcript-aware compatibility against non-fusion contig structure
- consider a read as contrary evidence if it supports a non-fusion interpretation across the relevant partner boundary

### 2. FAR values are still provisional for long reads

Because LR contrary-read support is incomplete, `FAR_left` and `FAR_right` should currently be treated as provisional for long-read runs.

### 3. No dedicated LR retrieval/visualization evidence path yet

The current implementation focuses on getting the correct long-read-supported fusion calls through the FI characterization pipeline.

It does not yet introduce a full LR-specific replacement for the old spanning-read retrieval machinery.


## Suggested Regression Checks

When revisiting this implementation on other datasets, verify:

1. `fusions.list` matches the intended candidate set from `ctat-LR-fusion`
2. final FI fusion names match all intended long-read candidates
3. no `EMadj` files are created for `--read_type long`
4. `est_J == JunctionReadCount` and `est_S == 0.00`
5. downstream files are still produced:
   - post blast/promiscuity filter
   - FFPM output
   - microhomology output
   - annotated final fusion report
6. inspect a few representative events for:
   - breakpoint preservation
   - splice-type consistency
   - counter-read plausibility


## Example Tested Command

```bash
/home/unix/bhaas/GITHUB/CTAT_FUSIONS/FusionInspector/FusionInspector \
  --left_fq fus_ev_reads.fq \
  --CPU 7 \
  --read_type long \
  --fusions fusions.list \
  --min_per_id 70 \
  -O FI-min70pid-lrcompat-noEM-counters
```


## Current Recommendation

Use the LR-compatible no-EM path as the baseline for future long-read troubleshooting.

Best current reference output:

- `/home/unix/bhaas/projects/Melanoma/FUSION-reads-only/FI-min70pid-lrcompat-noEM-counters`

If future discrepancies appear, compare:

- the FI candidate list in `fusions.list`
- the LR support table in `fi_workdir/finspector.fusion_preds.coalesced.summary`
- the final FI report in `finspector.FusionInspector.fusions.tsv`

