# CUDLL BAMs to FASTQ Converter

This plugin processes two CUDLL BAM files and creates a single merged FASTQ file containing primary reads.

## Features

- Extracts primary reads (excludes secondary and supplementary alignments) from two BAM files
- Reads from the supplemental BAM take priority
- Any read names in the main BAM that also appear in the supplemental BAM are excluded
- Outputs a compressed FASTQ.gz file

## Components

### Python Script (`Docker/cudll_bams_to_fastq.py`)

The main processing script that:
1. Reads the supplemental BAM and extracts all primary read names
2. Writes all primary reads from the supplemental BAM to the output FASTQ
3. Writes primary reads from the main BAM that don't appear in the supplemental BAM
4. Reports statistics on the number of reads processed

### Dockerfile (`Docker/Dockerfile`)

Extends the FusionInspector Docker image (`trinityctat/fusioninspector:2.11.0`) and adds:
- `pysam` Python library for BAM file processing
- The `cudll_bams_to_fastq.py` script

### WDL Workflow (`WDL/cudll_bams_to_fastq.wdl`)

Terra-compatible WDL workflow that:
- Takes two BAM files and a sample name as input
- Runs the conversion process in a Docker container
- Outputs the merged FASTQ file and a log file

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
    --output /data/merged.fastq.gz
```

### Directly (if pysam is installed)

```bash
./Docker/cudll_bams_to_fastq.py \
  --main-bam main.bam \
  --supp-bam supp.bam \
  --output merged.fastq.gz
```

## Running on Terra

1. Upload the WDL file to Terra
2. Configure inputs:
   - `cudll_main_bam`: Path to the main CUDLL BAM file
   - `cudll_supp_bam`: Path to the supplemental CUDLL BAM file
   - `sample_name`: Sample identifier for output naming
   - `docker` (optional): Docker image to use (default: `trinityctat/fusioninspector-cudll:latest`)
3. Launch the workflow

## Outputs

- **merged_fastq**: Compressed FASTQ file containing merged reads
- **log_file**: Processing log with statistics

## Requirements

- BAM files must be readable by pysam
- BAM files should contain sequence and quality information
- Sufficient disk space for intermediate processing and output files
