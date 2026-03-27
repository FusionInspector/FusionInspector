# CUDLL BAMs to CB/UMI Table Converter

This plugin processes one or two CUDLL BAM files and creates a tab-delimited table containing read names with their corresponding cell barcodes (CB) and UMIs. It also emits a barcode-level summary of unique UMI counts and a PDF knee plot of barcode rank versus unique UMI count.

## Features

- Extracts primary reads (excludes secondary and supplementary alignments) from one or two BAM files
- Reads from the supplemental BAM take priority when provided
- Any read names in the main BAM that also appear in the supplemental BAM are excluded
- Extracts cell barcode from CB tag (configurable)
- Extracts UMI from XM tag (configurable)
- Outputs a compressed tab-delimited file (.tsv.gz)
- Outputs a barcode-level table of unique UMI counts per cell barcode
- Outputs a PDF knee plot for barcode rank by unique UMI count

## Output Format

The output is a tab-delimited file with the following columns:

```
read_name    cell_barcode    UMI
```

Reads missing the CB or UMI tags will have "NA" in the corresponding column.

## Barcode UMI Summary

The plugin also writes a second tab-delimited file with one row per cell barcode:

```
cell_barcode    umi_count
```

`umi_count` is the number of distinct UMIs observed for that barcode across the retained reads. Rows where the barcode or UMI tag is missing are excluded from this summary and from the knee plot.

## Components

### Python Script (`Docker/cudll_bams_to_cb_umi.py`)

The main processing script that:
1. Optionally reads the supplemental BAM and extracts all primary read names
2. Writes all primary reads from the supplemental BAM with their CB/UMI tags to the output table when provided
3. Writes primary reads from the main BAM (excluding reads already present in the supplemental BAM, if provided) with their CB/UMI tags
4. Tracks distinct barcode/UMI pairs and writes a barcode summary table
5. Generates a PDF knee plot of barcode rank versus unique UMI count
6. Reports statistics on the number of reads processed and missing tags

### Dockerfile (`Docker/Dockerfile`)

Extends the FusionInspector Docker image (`trinityctat/fusioninspector:2.11.1`) and adds:
- `pysam` Python library for BAM file processing
- `matplotlib` for PDF knee-plot generation
- The `cudll_bams_to_cb_umi.py` script

### WDL Workflow (`WDL/cudll_bams_to_cb_umi.wdl`)

Terra-compatible WDL workflow that:
- Takes two BAM files and a sample name as input
- Allows configuration of CB and UMI tag names (defaults: CB and XM)
- Runs the conversion process in a Docker container
- Outputs the CB/UMI table, barcode summary table, knee-plot PDF, and a log file

## Building the Docker Image

```bash
cd Docker
./build_docker.sh
```

Or manually:

```bash
cd Docker
docker build -t trinityctat/cudll-to-cb-umi:latest .
```

## Pushing the Docker Image

```bash
cd Docker
./push_docker.sh
```

## Running Locally

### Using Docker

```bash
docker run --rm -v /path/to/data:/data \
  trinityctat/cudll-to-cb-umi:latest \
  cudll_bams_to_cb_umi.py \
    --main-bam /data/main.bam \
    --output /data/sample.cb_umi.tsv.gz \
    --barcode-umi-counts /data/sample.barcode_umi_counts.tsv \
    --knee-plot-pdf /data/sample.barcode_knee_plot.pdf \
    --cb-tag CB \
    --umi-tag XM
```

### Directly (if pysam is installed)

```bash
./Docker/cudll_bams_to_cb_umi.py \
  --main-bam main.bam \
  --output sample.cb_umi.tsv.gz \
  --barcode-umi-counts sample.barcode_umi_counts.tsv \
  --knee-plot-pdf sample.barcode_knee_plot.pdf \
  --cb-tag CB \
  --umi-tag XM
```

## Running on Terra

1. Upload the WDL file to Terra
2. Configure inputs:
   - `cudll_main_bam`: Path to the main CUDLL BAM file
   - `cudll_supp_bam` (optional): Path to the supplemental CUDLL BAM file
   - `sample_name`: Sample identifier for output naming
   - `cb_tag` (optional): BAM tag for cell barcode (default: CB)
   - `umi_tag` (optional): BAM tag for UMI (default: XM)
   - `docker` (optional): Docker image to use (default: `trinityctat/cudll-to-cb-umi:latest`)
3. Launch the workflow

## Execution Test

A deterministic local execution test is available under `testing/`.

```bash
cd plugins/CUDLL_bams_to_CB_UMI
./testing/run_execution_test.sh
```

This generates synthetic BAMs, runs the plugin, and validates the read-level table, barcode-level UMI counts, and knee-plot PDF.

## Outputs

- **cb_umi_table**: Compressed tab-delimited file containing read_name, cell_barcode, and UMI
- **barcode_umi_counts**: Tab-delimited file containing cell_barcode and distinct UMI count
- **barcode_knee_plot_pdf**: PDF plot of barcodes ranked by distinct UMI count
- **log_file**: Processing log with statistics and warnings for missing tags

## Requirements

- BAM files must be readable by pysam
- BAM files should contain the specified CB and UMI tags
- Sufficient disk space for intermediate processing and output files

## Parameters

- `--main-bam`: Main CUDLL BAM file (required)
- `--supp-bam`: Supplemental CUDLL BAM file that takes priority (optional)
- `--output`: Output file path, can be .gz compressed (required)
- `--barcode-umi-counts`: Output TSV for unique UMI counts per barcode (required)
- `--knee-plot-pdf`: Output PDF for barcode knee plot (required)
- `--cb-tag`: BAM tag for cell barcode (default: CB)
- `--umi-tag`: BAM tag for UMI (default: XM)
