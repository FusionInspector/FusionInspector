# Long Read Support for FusionInspector

## Overview

FusionInspector now supports PacBio and Oxford Nanopore long reads using minimap2 as an alternative aligner, while maintaining full backward compatibility with existing STAR/Illumina workflows.

## Implementation Summary

### New Features

1. **Minimap2 Integration**: Added minimap2 aligner support for efficient long read alignment
2. **Read Type Detection**: Automatic aligner selection based on read type
3. **Single-End Support**: Enhanced single-end read handling for long reads
4. **Backward Compatible**: All existing STAR/Illumina workflows remain unchanged

### Key Components Modified

#### 1. New Files Created
- **`util/run_FI_minimap2.pl`**: Minimap2 wrapper script (mirrors run_FI_STAR.pl architecture)
  - Builds minimap2 index (.mmi file)
  - Converts GTF to junction BED using paftools.js
  - Handles genome patching (concatenates reference + fusion contigs)
  - Supports fusion-contigs-only and full-genome modes
  - Manages read groups for samples files
  - Post-filters BAM for fusion-only reads when requested

#### 2. Main Script Updates (`FusionInspector`)
- Added three new command-line arguments:
  - `--aligner [STAR|minimap2]`: Select aligner (default: STAR)
  - `--read_type [short|long]`: Specify read type (default: short)
  - `--minimap2_params "..."`: Extra minimap2 parameters (default: "-x splice")
- Added validation logic:
  - Long reads require single-end input only (no --right_fq)
  - Automatically sets minimap2 for long reads
  - Warns if using minimap2 with short reads
- Updated alignment section to be aligner-agnostic:
  - Dynamic aligner script selection
  - Dynamic BAM file naming based on aligner
  - Aligner-specific parameter handling

#### 3. Downstream Compatibility Fixes
- **`util/get_fusion_JUNCTION_reads_from_fusion_contig_bam.pl`**:
  - Added NH:i: tag fallback (minimap2 with --secondary=no doesn't output this tag)
  - Defaults to NH:i:1 when tag is absent

- **`util/get_fusion_SPANNING_reads_from_bam.from_chim_summary.pl`**:
  - Added NH:i: tag fallback (same as above)

#### 4. Docker Updates
- **`Docker/Dockerfile`**:
  - Already had minimap2 v2.26 installed
  - Already had k8 JavaScript shell installed
  - Added paftools.js to $BIN for GTF to BED conversion

## Usage

### Basic Long Read Usage

#### PacBio HiFi Reads
```bash
FusionInspector \
    --fusions fusion_targets.txt \
    --genome_lib_dir $CTAT_GENOME_LIB \
    --left_fq pacbio_hifi.fastq.gz \
    --read_type long \
    --out_dir FI_hifi \
    --CPU 8
```

The `--read_type long` automatically selects minimap2 as the aligner.

#### PacBio HiFi with Explicit Aligner
```bash
FusionInspector \
    --fusions fusion_targets.txt \
    --genome_lib_dir $CTAT_GENOME_LIB \
    --left_fq pacbio_hifi.fastq.gz \
    --aligner minimap2 \
    --minimap2_params "-x splice:hq" \
    --out_dir FI_hifi \
    --CPU 8
```

Note: `-x splice:hq` is recommended for PacBio HiFi reads (high quality).

#### Oxford Nanopore Reads
```bash
FusionInspector \
    --fusions fusion_targets.txt \
    --genome_lib_dir $CTAT_GENOME_LIB \
    --left_fq ont_reads.fastq.gz \
    --read_type long \
    --out_dir FI_ont \
    --CPU 8
```

#### Fusion Contigs Only Mode (Faster)
```bash
FusionInspector \
    --fusions fusion_targets.txt \
    --genome_lib_dir $CTAT_GENOME_LIB \
    --left_fq pacbio_hifi.fastq.gz \
    --read_type long \
    --fusion_contigs_only \
    --out_dir FI_hifi_contigs \
    --CPU 8
```

### Traditional Short Read Usage (Unchanged)

#### Illumina Paired-End (Default)
```bash
FusionInspector \
    --fusions fusion_targets.txt \
    --genome_lib_dir $CTAT_GENOME_LIB \
    --left_fq reads_1.fastq.gz \
    --right_fq reads_2.fastq.gz \
    --out_dir FI \
    --CPU 8
```

No changes needed - STAR is still the default aligner for short reads.

## Technical Details

### Alignment Modes

FusionInspector supports two alignment modes (same for both STAR and minimap2):

1. **Default Mode (Recommended for Long Reads)**
   - Aligns to full reference genome + fusion contigs
   - Uses combined GTF (reference annotations + fusion contig annotations)
   - Better context for long read alignment
   - Reduces false positives
   - Command: FusionInspector with minimap2 (no `--fusion_contigs_only` flag)

2. **Fusion Contigs Only Mode**
   - Aligns to fusion contigs mini-genome only
   - Uses fusion contig GTF only
   - Faster, more focused on fusion-supporting reads
   - Command: FusionInspector with `--fusion_contigs_only` flag

### Minimap2 Parameters

Default minimap2 command:
```bash
minimap2 -ax splice -t $CPU -G $max_mate_dist \
  -uf --MD -L --junc-bed $splice_bed \
  $mm2_index $reads | \
  samtools sort -@ $CPU -m 4G -o $output.bam -
```

Key parameters:
- `-ax splice`: Splice-aware alignment for RNA reads
- **Secondary alignments allowed** (default): Critical for reads that align to both genome and fusion contigs
  - Non-fusion transcripts need to multi-map to genomic locations AND fusion contigs
  - The NH:i: tag tracks alignment count for downstream filtering
  - Use `-N <int>` to control max secondary alignments (default: 5)
- `-uf`: Align with forward strand direction
- `--MD`: Output MD tag for mismatches
- `-L`: Long CIGAR format
- `--junc-bed`: Splice junction hints from GTF

For PacBio HiFi, consider using `-x splice:hq` for higher quality settings.

### Genome Patching

Unlike STAR which supports on-the-fly patching via `--genomeFastaFiles`, minimap2 requires:
1. Concatenating genome + fusion contigs: `cat $genome $patch > $combined.fa`
2. Concatenating GTFs: `cat $genome.gtf $patch.gtf > $combined.gtf`
3. Building index on combined genome: `minimap2 -d $combined.mmi $combined.fa`

This is handled automatically by `run_FI_minimap2.pl`.

### Read Group Handling

For samples files, minimap2 doesn't have direct read group support like STAR's `--outSAMattrRGline`.
The script uses `samtools addreplacerg` to add read groups post-alignment.

### Fusion Read Filtering

For `--only_fusion_reads` mode, the script post-filters the BAM using AWK:
```bash
samtools view -h $bam | \
  awk 'BEGIN{OFS="\t"} /^@/ {print; next} $3 ~ /--/ {print}' | \
  samtools view -bo $filtered.bam -
```

This filters for reference names containing "--" (fusion contig naming convention).

## Compatibility Notes

### SAM/BAM Format
- Both STAR and minimap2 output standard SAM/BAM format
- Existing SAM_reader.pm and SAM_entry.pm work with both aligners
- Single-end read handling already supported throughout the pipeline

### NH:i: Tag
- Both STAR and minimap2 output NH:i:X tag indicating number of alignments
- This tag is used by downstream scripts to assess alignment uniqueness
- **Important**: Secondary alignments must be allowed for proper multi-mapping behavior
- Downstream scripts include a safety fallback to NH:i:1 if tag is absent

### Error Profiles
Long reads have different error profiles than Illumina:
- PacBio CLR: ~15% error rate
- PacBio HiFi: ~0.1% error rate (similar to Illumina)
- ONT: ~5-15% error rate (improving with newer basecallers)

The `--min_per_id` parameter (default: 96) may need adjustment for long reads:
```bash
--min_per_id 90  # More permissive for PacBio CLR or ONT
--min_per_id 96  # Default, suitable for PacBio HiFi
```

## Testing

### Unit Testing
Test minimap2 wrapper independently:
```bash
cd /path/to/FusionInspector/test
util/run_FI_minimap2.pl \
  --genome data/minigenome.fa \
  --reads test.reads_1.fastq.gz \
  -G data/minigenome.gtf \
  --CPU 2 \
  --out_prefix test_mm2
```

### Integration Testing
Run with test data:
```bash
cd /path/to/FusionInspector/test
./FusionInspector \
  --fusions fusion_targets.A.txt \
  --genome_lib_dir $CTAT_GENOME_LIB \
  --left_fq test.readsA_1.fastq.gz \
  --aligner minimap2 \
  --out_dir test_mm2_integration
```

### Backward Compatibility Testing
Ensure existing tests still pass:
```bash
cd /path/to/FusionInspector/test
make test  # Should use STAR by default
```

## Performance Expectations

Minimap2 is typically:
- **Faster** than STAR for long reads (5-10x speedup)
- **Lower memory** usage (no large genome index directory)
- **Smaller index** (.mmi file vs. STAR's multi-file index)

Example timings (approximate):
- STAR with Illumina PE: ~10-20 minutes for 10M read pairs
- Minimap2 with PacBio: ~5-10 minutes for 100K long reads

## Troubleshooting

### paftools.js Not Found
If you see warnings about paftools.js:
```bash
Warning: paftools.js not found. Skipping junction BED creation.
```

This is usually not critical, but for best results:
1. Ensure k8 JavaScript shell is installed: `which k8`
2. Ensure paftools.js is in PATH: `which paftools.js`
3. Or update Docker image to include paftools.js (already done in latest Dockerfile)

### Long Reads with --right_fq Error
```
Error, long reads (--read_type long) are single-end only. Cannot specify --right_fq.
```

Long reads are sequenced as single molecules, not paired-end. Remove `--right_fq` argument.

### Minimap2 Not in PATH
```
Error, cannot locate minimap2 program. Be sure it's in your PATH setting.
```

Install minimap2:
```bash
# Via conda
conda install -c bioconda minimap2

# Or download binary
curl -L https://github.com/lh3/minimap2/releases/download/v2.26/minimap2-2.26_x64-linux.tar.bz2 | tar -jxvf -
cp minimap2-2.26_x64-linux/minimap2 /usr/local/bin/
```

## Future Enhancements

Potential future improvements:
1. Automatic read type detection from FASTQ (based on read length)
2. Long read-specific quality metrics
3. Integration with long-read fusion detection tools (e.g., CTAT-LR-fusion)
4. Support for additional long-read aligners (e.g., pbmm2)

## References

- Minimap2: https://github.com/lh3/minimap2
- CTAT-LR-fusion: https://github.com/NCIP/CTAT-LR-fusion
- FusionInspector: https://github.com/FusionInspector/FusionInspector

## Version History

- **v2.x.x** (2026-03): Added long read support with minimap2
  - New `--aligner`, `--read_type`, and `--minimap2_params` arguments
  - Created `run_FI_minimap2.pl` wrapper script
  - Enhanced NH:i: tag compatibility for downstream scripts
  - Updated Docker to include paftools.js
