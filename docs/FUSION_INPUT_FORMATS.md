# FusionInspector: Supported Fusion Input Formats

## Overview

FusionInspector now supports multiple input formats for the `--fusions` parameter:

1. **Simple list format** (original)
2. **STAR-Fusion output format** (new)
3. **CTAT-LR-Fusion output format** (new)

The tool automatically detects the input format and extracts the fusion targets accordingly.

## Format 1: Simple List (Original)

A plain text file with one fusion per line in `geneA--geneB` format:

```
# Optional comment lines starting with #
ACACA--STAC2
EML4--ALK
BCR--ABL1
```

**Format:**
- One fusion per line
- Fusion format: `geneA--geneB` or `geneA::geneB`
- Lines starting with `#` are treated as comments
- Empty lines are ignored

## Format 2: STAR-Fusion Output

The full output file from STAR-Fusion with the `#FusionName` column:

```
#FusionName	JunctionReadCount	SpanningFragCount	...
EML4--ALK	150	50	...
ACACA--STAC2	100	30	...
EML4--ALK	80	20	...  (different breakpoint)
BCR--ABL1	200	60	...
```

**How it works:**
- FusionInspector detects the `#FusionName` header column
- Extracts fusion names from the `#FusionName` column
- **Automatically deduplicates** fusions (same fusion with different breakpoints → single entry)
- Creates a simple list file for downstream processing

**Example:**
```bash
FusionInspector \
  --fusions star-fusion.fusion_predictions.tsv \
  --genome_lib_dir $CTAT_GENOME_LIB \
  --left_fq reads_1.fq.gz \
  --right_fq reads_2.fq.gz \
  --out_dir FI_output
```

## Format 3: CTAT-LR-Fusion Output

The full output file from CTAT-LR-Fusion with the `#FusionName` column:

```
#FusionName	num_LR	LR_FFPM	...
EML4--ALK	25	5.2	...
ACACA--STAC2	15	3.1	...
BCR--ABL1	30	6.5	...
```

**How it works:**
- Same as STAR-Fusion (detects `#FusionName` column)
- Extracts and deduplicates fusion names
- Processes identically regardless of short-read or long-read origin

**Example:**
```bash
FusionInspector \
  --fusions ctat-LR-fusion.fusions.tsv \
  --genome_lib_dir $CTAT_GENOME_LIB \
  --left_fq long_reads.fq.gz \
  --read_type long \
  --out_dir FI_output
```

## Multiple Input Files

You can provide multiple fusion files (comma-delimited, no spaces):

```bash
FusionInspector \
  --fusions star_fusion.tsv,ctat_lr_fusion.tsv,manual_list.txt \
  --genome_lib_dir $CTAT_GENOME_LIB \
  --left_fq reads.fq.gz \
  --out_dir FI_output
```

**Behavior:**
- Each file is processed independently
- Files can be in different formats (mixed)
- All extracted fusions are combined
- Duplicates across files are handled automatically

## Automatic Format Detection

FusionInspector automatically detects the format by:

1. Reading the first line of the file
2. If it starts with `#` and contains `#FusionName` → **STAR-Fusion/CTAT-LR-Fusion format**
3. Otherwise → **Simple list format**

No user action required!

## Deduplication Logic

When using STAR-Fusion or CTAT-LR-Fusion output:

### Why Deduplication?

The same fusion can appear multiple times with:
- Different breakpoints
- Different isoforms
- Different evidence types

**Example:**
```
#FusionName	LeftBreakpoint	RightBreakpoint
EML4--ALK	chr2:42522656	chr2:29415640
EML4--ALK	chr2:42522700	chr2:29415640
EML4--ALK	chr2:42525234	chr2:29415640
```

All three represent the same `EML4--ALK` fusion.

### What Happens

FusionInspector extracts unique fusion names:
- Input: 3 lines for `EML4--ALK`
- Output: 1 entry for `EML4--ALK`

FusionInspector will then:
1. Build fusion contigs for `EML4--ALK` (once)
2. Align reads to the fusion contig
3. Identify **all** breakpoints during alignment (including the 3 from above)

## Preprocessing Details

### Temporary Files

When STAR-Fusion/CTAT-LR-Fusion format is detected:

1. A simple list file is created: `<output_dir>/fi_workdir/fusion_targets.txt`
2. Contains unique fusion names (deduplicated)
3. Used internally for downstream processing

**Example content:**
```
ACACA--STAC2
BCR--ABL1
EML4--ALK
```

### Log Messages

When processing, you'll see:

```
INFO: Detected STAR-Fusion/CTAT-LR-Fusion output format in star-fusion.tsv
INFO: Extracted 3 unique fusions from star-fusion.tsv
INFO: Wrote simple fusion list to <output_dir>/fi_workdir/fusion_targets.txt
```

Or for simple list:

```
INFO: Using fusion file fusions.txt as-is (simple list format)
```

## Gzip Support

All formats support gzip compression:

```bash
# All of these work
--fusions fusions.txt
--fusions fusions.txt.gz
--fusions star-fusion.fusion_predictions.tsv.gz
--fusions ctat-lr-fusion.fusions.tsv.gz
```

## Error Handling

### No Fusions Found

If no valid fusions are detected:

```
Warning: No fusions found in #FusionName column of <file>
No valid fusions found after preprocessing input files.
```

The program exits gracefully.

### Missing #FusionName Column

If file has header but no `#FusionName` column:

```
Error: Could not find #FusionName column in <file>
```

## Best Practices

### For STAR-Fusion Users

Use the fusion predictions file directly:

```bash
FusionInspector \
  --fusions star-fusion.fusion_predictions.tsv \
  ...
```

**Advantages:**
- No manual extraction needed
- Automatic deduplication
- Preserves all fusion candidates

### For CTAT-LR-Fusion Users

Use the fusions output file directly:

```bash
FusionInspector \
  --fusions ctat-LR-fusion.fusions.tsv \
  --read_type long \
  ...
```

### For Custom Lists

Create a simple text file:

```bash
# fusions_to_validate.txt
GENE1--GENE2
GENE3--GENE4
```

```bash
FusionInspector \
  --fusions fusions_to_validate.txt \
  ...
```

## Compatibility

- **Backward compatible**: Existing simple list workflows unchanged
- **Forward compatible**: New format detection is automatic
- **Version**: Available in v2.11.0+

## Summary

| Input Format | Detection | Deduplication | Example File |
|--------------|-----------|---------------|--------------|
| Simple list | No `#FusionName` header | N/A | `fusions.txt` |
| STAR-Fusion | `#FusionName` column | Yes, by fusion name | `star-fusion.fusion_predictions.tsv` |
| CTAT-LR-Fusion | `#FusionName` column | Yes, by fusion name | `ctat-LR-fusion.fusions.tsv` |

**Key takeaway:** Just pass your fusion file directly - FusionInspector handles the rest!
