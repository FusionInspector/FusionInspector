# CUDLL BAMs to FASTQ Converter

This plugin processes two CUDLL BAM files and creates a single merged FASTQ file containing primary reads.

## Features

- Extracts primary reads (excludes secondary and supplementary alignments) from two BAM files
- Reads from the supplemental BAM take priority
- Any read names in the main BAM that also appear in the supplemental BAM are excluded
- Formats read names as: `cell_barcode^umi^read_name`
- Extracts cell barcode from CB tag (configurable)
- Extracts UMI from XM tag (configurable)
- Outputs a compressed FASTQ.gz file

## Components

### Python Script (`Docker/cudll_bams_to_fastq.py`)

The main processing script that:
1. Reads the supplemental BAM and extracts all primary read names
2. Writes all primary reads from the supplemental BAM to the output FASTQ with formatted names
3. Writes primary reads from the main BAM that don't appear in the supplemental BAM
4. Formats read names as `cell_barcode^umi^read_name` using CB and UMI tags
5. Reports statistics on the number of reads processed and missing tags

### Dockerfile (`Docker/Dockerfile`)

Extends the FusionInspector Docker image (`trinityctat/fusioninspector:2.11.1`) and adds:
- `pysam` Python library for BAM file processing
- The `cudll_bams_to_fastq.py` script

### WDL Workflow (`WDL/cudll_bams_to_fastq.wdl`)

Terra-compatible WDL workflow that:
- Takes two BAM files and a sample name as input
- Allows configuration of CB and UMI tag names (defaults: CB and XM)
- Runs the conversion process in a Docker container
- Outputs the merged FASTQ file with formatted read names and a log file

## Building the Docker Image

```bash
cd Docker
docker build -t trinityctat/fusioninspector-cudll:latest .
```

## Running Locally

### Using Docker

```bash
docker run --rm -v /path/to/data:/data \
  trinityctat/fusioninspector-cudll:latest \
  cudll_bams_to_fastq.py \
    --main-bam /data/main.bam \
    --supp-bam /data/supp.bam \
    --output /data/merged.fastq.gz \
    --cb-tag CB \
    --umi-tag XM
```

### Directly (if pysam is installed)

```bash
./Docker/cudll_bams_to_fastq.py \
  --main-bam main.bam \
  --supp-bam supp.bam \
  --output merged.fastq.gz \
  --cb-tag CB \
  --umi-tag XM
```

## Running on Terra

1. Upload the WDL file to Terra
2. Configure inputs:
   - `cudll_main_bam`: Path to the main CUDLL BAM file
   - `cudll_supp_bam`: Path to the supplemental CUDLL BAM file
   - `sample_name`: Sample identifier for output naming
   - `cb_tag` (optional): BAM tag for cell barcode (default: CB)
   - `umi_tag` (optional): BAM tag for UMI (default: XM)
   - `docker` (optional): Docker image to use (default: `trinityctat/fusioninspector-cudll:latest`)
3. Launch the workflow

## Outputs

- **merged_fastq**: Compressed FASTQ file containing merged reads with formatted names (`cell_barcode^umi^read_name`)
- **log_file**: Processing log with statistics and warnings for missing tags

## Read Name Format

FASTQ read names are formatted as: `cell_barcode^umi^read_name`

For example:
```
@AAACCCAAGAAACACT-1^TTTGGGAAA^A00123:123:H5YM3DSXY:1:1101:1234:1000/1
```

Reads missing CB or UMI tags will have "NA" in the corresponding position.

## Requirements

- BAM files must be readable by pysam
- BAM files should contain sequence and quality information
- BAM files should contain CB and UMI tags (or alternative tags as specified)
- Sufficient disk space for intermediate processing and output files

## Parameters

- `--main-bam`: Main CUDLL BAM file (required)
- `--supp-bam`: Supplemental CUDLL BAM file that takes priority (required)
- `--output`: Output FASTQ.gz file path (required)
- `--cb-tag`: BAM tag for cell barcode (default: CB)
- `--umi-tag`: BAM tag for UMI (default: XM)
