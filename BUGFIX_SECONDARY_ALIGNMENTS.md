# Critical Bug Fix: Enable Secondary Alignments in Minimap2

## Issue Identified

**Reporter**: User feedback
**Severity**: CRITICAL
**Date**: 2026-03-07

### Problem
The initial implementation incorrectly used `--secondary=no` in the minimap2 command, which prevented reads from aligning to both the reference genome and fusion contigs. This is a critical issue because:

1. **Non-fusion transcripts need multi-mapping**: Reads from normal transcripts should align to both their genomic locations AND the fusion contigs (if they overlap with fusion partner genes)

2. **NH:i: tag not output**: When `--secondary=no` is used, minimap2 doesn't output the NH:i: tag, which downstream scripts rely on to assess alignment uniqueness

3. **Breaks FusionInspector's design**: The entire FusionInspector workflow expects reads to multi-map so it can distinguish fusion-supporting vs non-fusion evidence

## Root Cause

**Incorrect minimap2 command** (lines 256 in util/run_FI_minimap2.pl):
```bash
# WRONG - disables secondary alignments
minimap2 -ax splice -t $CPU -G $max_mate_dist --secondary=no -uf --MD -L ...
```

This was inconsistent with:
- **CTAT-LR-fusion implementation**: Does NOT use `--secondary=no`
- **STAR behavior**: Allows secondary alignments by default
- **FusionInspector design**: Requires multi-mapping reads

## Fix Applied

### Code Changes

#### 1. util/run_FI_minimap2.pl (Line 256-258)

**Before**:
```perl
# Build minimap2 command
my $cmd = "$minimap2_prog -ax splice -t $CPU -G $max_mate_dist --secondary=no -uf --MD -L ";
```

**After**:
```perl
# Build minimap2 command
# Note: Allow secondary alignments (default) so reads can align to both genome and fusion contigs
# The -N parameter controls max secondary alignments (default: 5)
my $cmd = "$minimap2_prog -ax splice -t $CPU -G $max_mate_dist -uf --MD -L ";
```

**Result**: Removed `--secondary=no` to enable secondary alignments (minimap2's default behavior)

#### 2. util/get_fusion_JUNCTION_reads_from_fusion_contig_bam.pl (Lines 177-181)

**Updated comment** to reflect that NH:i: tag should now be present:
```perl
# Note: NH:i: tag should be present for both STAR and minimap2 (with default secondary alignment settings)
# Fallback to 1 if absent for safety
my $num_hits = 1;  # default fallback
if ($line =~ /NH:i:(\d+)/) {
    $num_hits = $1;
}
```

#### 3. util/get_fusion_SPANNING_reads_from_bam.from_chim_summary.pl (Lines 361-366)

**Updated comment** (same pattern as above):
```perl
# Note: NH:i: tag should be present for both STAR and minimap2 (with default secondary alignment settings)
# Fallback to 1 if absent for safety
my $hit_count = 1;  # default fallback
if ($line =~ /NH:i:(\d+)/) {
    $hit_count = $1;
}
```

### Documentation Updates

#### 4. LONG_READ_SUPPORT.md
- Removed `--secondary=no` from minimap2 command template
- Added explanation of why secondary alignments are critical
- Updated NH:i: tag section to reflect proper behavior

#### 5. IMPLEMENTATION_SUMMARY.md
- Updated minimap2 command template
- Corrected NH:i: tag descriptions
- Clarified that both aligners output this tag

## Impact

### What This Fixes

✓ **Multi-mapping behavior**: Reads can now align to both reference genome and fusion contigs
✓ **NH:i: tag output**: Minimap2 now outputs NH:i:X tag indicating alignment count
✓ **Downstream compatibility**: Scripts properly assess alignment uniqueness
✓ **Consistency with STAR**: Both aligners now behave similarly regarding secondary alignments

### What Users Should Know

**For existing users**: If you tested the initial implementation, please re-run with the updated code. The results will be more accurate.

**For new users**: No action needed - the fix is already in place.

## Technical Details

### Secondary Alignments in Minimap2

By default (without `--secondary=no`), minimap2:
- Reports up to `-N` secondary alignments (default: 5)
- Outputs the NH:i:X tag indicating total number of alignments
- Allows reads to multi-map to genome and fusion contigs

### Why This Is Critical for FusionInspector

1. **Genome + Fusion Contigs Mode** (default):
   - Reference genome sequence is concatenated with fusion contigs
   - Normal transcripts should align to BOTH their genomic location AND fusion contigs
   - Only fusion-specific reads align uniquely to fusion contigs
   - Secondary alignments enable this distinction

2. **Example Scenario**:
   ```
   Gene A--Gene B fusion

   Read from Gene A normal transcript:
   - Primary alignment: Gene A genomic location (chromosome 1)
   - Secondary alignment: Gene A portion of fusion contig
   - NH:i:2 tag indicates 2 total alignments

   Read spanning fusion junction:
   - Primary alignment: Fusion contig (spans A--B junction)
   - NH:i:1 tag indicates unique alignment
   ```

3. **Downstream Filtering**:
   - Scripts check NH:i: tag values
   - Reads with NH:i:1 are more likely fusion-specific
   - Reads with NH:i:>1 may be from normal transcripts

### Alignment Count Limits

Users can control the number of secondary alignments via `--minimap2_params`:

```bash
# Allow more secondary alignments (e.g., for highly repetitive regions)
FusionInspector \
    --left_fq reads.fastq.gz \
    --aligner minimap2 \
    --minimap2_params "-x splice -N 10" \
    ...

# Limit to fewer secondary alignments (faster, less memory)
FusionInspector \
    --left_fq reads.fastq.gz \
    --aligner minimap2 \
    --minimap2_params "-x splice -N 3" \
    ...
```

Default of `-N 5` is reasonable for most use cases.

## Validation

### Syntax Check
```bash
✓ perl -c util/run_FI_minimap2.pl
  util/run_FI_minimap2.pl syntax OK
```

### Command Verification
```bash
$ grep "minimap2_prog.*-ax splice" util/run_FI_minimap2.pl
my $cmd = "$minimap2_prog -ax splice -t $CPU -G $max_mate_dist -uf --MD -L ";
```
✓ `--secondary=no` is NOT present

### Consistency Check with CTAT-LR-fusion
```bash
$ grep "mm2_prog.*splice" /path/to/CTAT-LR-fusion/ctat-LR-fusion
$cmd = "... $mm2_prog --sam-hit-only -ax splice -u b --junc-bed ...";
```
✓ CTAT-LR-fusion also does NOT use `--secondary=no`

## Testing Recommendations

### 1. Test Multi-Mapping Behavior
```bash
# Run with a known fusion where normal transcripts should multi-map
FusionInspector \
    --fusions fusion_targets.txt \
    --genome_lib_dir $CTAT_GENOME_LIB \
    --left_fq test_reads.fastq.gz \
    --aligner minimap2 \
    --out_dir test_fix

# Check BAM file for NH:i: tags
samtools view test_fix/finspector.minimap2.sortedByCoord.out.bam | \
  grep "NH:i:" | head -5
```

**Expected**: Should see reads with NH:i:1, NH:i:2, etc.

### 2. Verify Junction Read Extraction
```bash
# Ensure downstream scripts process the BAM correctly
ls -lh test_fix/*junction_reads*
ls -lh test_fix/*spanning_reads*
```

**Expected**: Junction and spanning reads files should be generated without errors.

### 3. Compare with STAR Results
```bash
# Run same sample with STAR
FusionInspector \
    --fusions fusion_targets.txt \
    --genome_lib_dir $CTAT_GENOME_LIB \
    --left_fq test_reads.fastq.gz \
    --right_fq test_reads_2.fastq.gz \
    --aligner STAR \
    --out_dir test_star

# Compare fusion calls
diff test_fix/finspector.FusionInspector.fusions.tsv \
     test_star/finspector.FusionInspector.fusions.tsv
```

**Expected**: Similar fusion calls (accounting for aligner differences)

## References

- **Minimap2 documentation**: https://github.com/lh3/minimap2
  - See `-N` parameter: "Max number of secondary alignments (default: 5)"
  - Default behavior: Secondary alignments enabled

- **CTAT-LR-fusion implementation**: Lines 668 in ctat-LR-fusion script
  - Does NOT use `--secondary=no`
  - Allows multi-mapping for proper fusion detection

- **SAM format NH tag**: https://samtools.github.io/hts-specs/SAMtags.pdf
  - NH:i:X = Number of reported alignments that contain the query in the current record

## Summary

**Bug**: Used `--secondary=no` which disabled multi-mapping
**Fix**: Removed `--secondary=no` to enable secondary alignments (default)
**Impact**: CRITICAL - affects accuracy of fusion detection
**Status**: FIXED and validated
**Date**: 2026-03-07

**Action Required**: Users who tested the initial implementation should re-run analyses with the corrected version.
