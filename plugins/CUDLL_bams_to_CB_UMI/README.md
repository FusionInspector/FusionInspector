# CUDLL BAMs to CB/UMI Table Converter

This plugin processes one or two CUDLL BAM files and creates a tab-delimited table containing read names with their corresponding cell barcodes (CB) and UMIs. It can also emit a barcode-level summary of unique UMI counts and a PDF knee plot of barcode rank versus unique UMI count.

## Features

- Extracts primary reads (excludes secondary and supplementary alignments) from one or two BAM files
- When a supplemental BAM is provided, primary reads from both BAMs are merged and deduplicated by read name
- Extracts cell barcode from CB tag (configurable)
- Extracts UMI from XM tag (configurable)
- Outputs a compressed tab-delimited file (.tsv.gz)
- Optionally outputs a barcode-level table of unique UMI counts per cell barcode
- Optionally outputs a PDF knee plot for barcode rank by unique UMI count

## Output Format

The output is a tab-delimited file with the following columns:

```
read_name    cell_barcode    UMI
```

Reads missing the CB or UMI tags will have "NA" in the corresponding column.
When a supplemental BAM is provided, the final table is emitted in sorted `read_name` order after deduplication.

## Barcode UMI Summary

When requested, the plugin writes a second tab-delimited file with one row per cell barcode:

```
cell_barcode    umi_count
```

`umi_count` is the number of distinct UMIs observed for that barcode across the retained reads. Rows where the barcode or UMI tag is missing are excluded from this summary and from the knee plot.

When only the main BAM is provided, reads are streamed directly to the output with no sorting step — read names are unique within a single BAM so no deduplication is needed.

When a supplemental BAM is also provided, each BAM is written to its own temporary TSV and sorted independently, then the two sorted files are merged with `sort -m` (streaming merge, negligible memory) before deduplication. Sorting each BAM separately halves the data seen by each `sort` invocation compared to combining them first.

`--sort-threads` maps to `sort --parallel`, and `--sort-memory-per-thread` is multiplied by the thread count to derive the total `sort -S` memory budget for each sort invocation. If `sort` is killed by `SIGKILL` under memory pressure, the script retries automatically with progressively smaller per-thread memory settings down to `64M`.

## Components

### Python Script (`Docker/cudll_bams_to_cb_umi.py`)

The main processing script that:
1. When only the main BAM is provided: streams primary reads directly to the output (no sort needed)
2. When a supplemental BAM is also provided: writes each BAM to its own temporary TSV, sorts them independently, then merges with `sort -m` and deduplicates by read name
3. Optionally writes valid barcode/UMI pairs to a temporary TSV, externally sorts them, and streams barcode-level unique UMI counts
4. Optionally generates a PDF knee plot of barcode rank versus unique UMI count
5. Reports statistics on the number of reads processed, deduplicated, and missing tags

### Dockerfile (`Docker/Dockerfile`)

Extends the FusionInspector Docker image (`trinityctat/fusioninspector:2.11.1`) and adds:
- `pysam` Python library for BAM file processing
- `matplotlib` for PDF knee-plot generation
- GNU `sort` for scalable barcode/UMI aggregation
- The `cudll_bams_to_cb_umi.py` script

### WDL Workflow (`WDL/cudll_bams_to_cb_umi.wdl`)

Terra-compatible WDL workflow that:
- Takes two BAM files and a sample name as input
- Allows configuration of CB and UMI tag names (defaults: CB and XM)
- Can optionally skip barcode-summary and knee-plot generation
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
    --barcode-umi-counts /data/sample.barcode_umi_counts.tsv.gz \
    --knee-plot-pdf /data/sample.barcode_knee_plot.pdf \
    --max-knee-plot-points 5000 \
    --sort-threads 4 \
    --sort-memory-per-thread 2G \
    --cb-tag CB \
    --umi-tag XM
```

### Directly (if pysam is installed)

```bash
./Docker/cudll_bams_to_cb_umi.py \
  --main-bam main.bam \
  --output sample.cb_umi.tsv.gz \
  --barcode-umi-counts sample.barcode_umi_counts.tsv.gz \
  --knee-plot-pdf sample.barcode_knee_plot.pdf \
  --max-knee-plot-points 5000 \
  --sort-threads 4 \
  --sort-memory-per-thread 2G \
  --cb-tag CB \
  --umi-tag XM
```

To generate only the primary read-level CB/UMI table, omit `--barcode-umi-counts` and `--knee-plot-pdf`.

```bash
./Docker/cudll_bams_to_cb_umi.py \
  --main-bam main.bam \
  --output sample.cb_umi.tsv.gz \
  --cb-tag CB \
  --umi-tag XM
```

## Running on Terra

1. Upload the WDL file to Terra
2. Configure inputs:
   - `cudll_main_bam`: Path to the main CUDLL BAM file
   - `cudll_supp_bam` (optional): Path to the supplemental CUDLL BAM file to merge before deduplicating by read name
   - `sample_name`: Sample identifier for output naming
   - `cb_tag` (optional): BAM tag for cell barcode (default: CB)
   - `umi_tag` (optional): BAM tag for UMI (default: XM)
   - `generate_summary_outputs` (optional): Set to `false` to emit only the primary CB/UMI table and log file
   - `sort_threads` (optional): Number of threads to use for each external sort (default: `2`)
   - `sort_memory_per_thread` (optional): Approximate RAM budget per sort thread, such as `512M` or `2G` (default: `2G`). If GNU `sort` is killed on large inputs, lowering this value can reduce peak memory pressure at the cost of more disk I/O.
   - `min_disk_gb` (optional): Minimum disk request for the task (default: `100`)
   - `extra_disk_gb` (optional): Fixed disk headroom added on top of the input-size-based estimate (default: `50`)
   - `disk_scale_factor` (optional): Multiplier applied to the combined input BAM size in GB when estimating disk (default: `6.0`)
   - `docker` (optional): Docker image to use (default: `trinityctat/cudll-to-cb-umi:latest`)
3. Launch the workflow

The WDL estimates task disk as:

```text
requested_disk_gb = max(
  min_disk_gb,
  ceil((size(main_bam_gb) + size(supp_bam_gb_if_present)) * disk_scale_factor) + extra_disk_gb
)
```

The multiplier is intentionally conservative because the workflow writes large temporary files for read-level sorting/deduplication and barcode/UMI aggregation.

The Terra runtime defaults target a `32 GB` machine with `cpu=8`, `sort_threads=2`, and `sort_memory_per_thread=2G`.

## Execution Test

A deterministic local execution test is available under `testing/`.

```bash
cd plugins/CUDLL_bams_to_CB_UMI
./testing/run_execution_test.sh
```

This generates synthetic BAMs, runs the plugin in both full-output and CB/UMI-only modes, and validates the outputs.

## Outputs

- **cb_umi_table**: Gzipped tab-delimited file containing read_name, cell_barcode, and UMI
- **barcode_umi_counts**: Optional gzipped tab-delimited file containing cell_barcode and distinct UMI count
- **barcode_knee_plot_pdf**: Optional PDF plot of barcodes ranked by distinct UMI count
- **log_file**: Processing log with statistics and warnings for missing tags

## Requirements

- BAM files must be readable by pysam
- BAM files should contain the specified CB and UMI tags
- Sufficient disk space for intermediate processing and output files

## Parameters

- `--main-bam`: Main CUDLL BAM file (required)
- `--supp-bam`: Supplemental CUDLL BAM file to merge before deduplicating by read name (optional)
- `--output`: Output file path, can be .gz compressed (required)
- `--barcode-umi-counts`: Output TSV for unique UMI counts per barcode (optional, but required if `--knee-plot-pdf` is provided)
- `--knee-plot-pdf`: Output PDF for barcode knee plot (optional, but required if `--barcode-umi-counts` is provided)
- `--max-knee-plot-points`: Maximum number of barcode ranks to render in the knee plot (default: 5000)
- `--sort-threads`: Number of threads to use for each external sort invocation (default: `2`)
- `--sort-memory-per-thread`: Approximate RAM budget per sort thread, for example `512M` or `2G` (default: `2G`). The script automatically retries with smaller budgets after a `SIGKILL`, but setting a lower value up front can still be useful on constrained Terra runtimes.
- `--cb-tag`: BAM tag for cell barcode (default: CB)
- `--umi-tag`: BAM tag for UMI (default: XM)
