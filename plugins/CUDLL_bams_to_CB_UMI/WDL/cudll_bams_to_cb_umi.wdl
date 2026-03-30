version 1.0

workflow CUDLLBamsToCbUmi {
    input {
        File cudll_main_bam
        File? cudll_supp_bam
        String sample_name
        String cb_tag = "CB"
        String umi_tag = "XM"
        Int max_knee_plot_points = 5000
        String docker = "trinityctat/cudll-to-cb-umi"
        Int cpu = 4
        Int memory_gb = 16
        Int disk_gb = 100
    }

    call BamsToCbUmi {
        input:
            main_bam = cudll_main_bam,
            supp_bam = cudll_supp_bam,
            sample_name = sample_name,
            cb_tag = cb_tag,
            umi_tag = umi_tag,
            max_knee_plot_points = max_knee_plot_points,
            docker = docker,
            cpu = cpu,
            memory_gb = memory_gb,
            disk_gb = disk_gb
    }

    output {
        File cb_umi_table = BamsToCbUmi.cb_umi_table
        File barcode_umi_counts = BamsToCbUmi.barcode_umi_counts
        File barcode_knee_plot_pdf = BamsToCbUmi.barcode_knee_plot_pdf
        File log_file = BamsToCbUmi.log_file
    }
}

task BamsToCbUmi {
    input {
        File main_bam
        File? supp_bam
        String sample_name
        String cb_tag
        String umi_tag
        Int max_knee_plot_points
        String docker
        Int cpu
        Int memory_gb
        Int disk_gb
    }

    String output_table = "~{sample_name}.cb_umi.tsv.gz"
    String barcode_umi_counts_file = "~{sample_name}.barcode_umi_counts.tsv"
    String knee_plot_pdf_file = "~{sample_name}.barcode_knee_plot.pdf"
    String log_file_name = "~{sample_name}.cudll_cb_umi.log"

    command <<<
        set -e
        set -o pipefail

        echo "Starting CUDLL BAM to CB/UMI table conversion" | tee ~{log_file_name}
        echo "Sample: ~{sample_name}" | tee -a ~{log_file_name}
        echo "Main BAM: ~{main_bam}" | tee -a ~{log_file_name}
        echo "Supp BAM: ~{if defined(supp_bam) then supp_bam else "none"}" | tee -a ~{log_file_name}
        echo "Output: ~{output_table}" | tee -a ~{log_file_name}
        echo "Barcode UMI counts: ~{barcode_umi_counts_file}" | tee -a ~{log_file_name}
        echo "Knee plot PDF: ~{knee_plot_pdf_file}" | tee -a ~{log_file_name}
        echo "CB Tag: ~{cb_tag}" | tee -a ~{log_file_name}
        echo "UMI Tag: ~{umi_tag}" | tee -a ~{log_file_name}
        echo "Max knee plot points: ~{max_knee_plot_points}" | tee -a ~{log_file_name}
        echo "---" | tee -a ~{log_file_name}

        script_help="$(cudll_bams_to_cb_umi.py --help 2>&1 || true)"

        if printf '%s\n' "${script_help}" | grep -q -- '--barcode-umi-counts'; then
            echo "Detected newer cudll_bams_to_cb_umi.py CLI; using native summary outputs" | tee -a ~{log_file_name}
            cudll_bams_to_cb_umi.py \
                --main-bam ~{main_bam} \
                ~{if defined(supp_bam) then "--supp-bam " + supp_bam else ""} \
                --output ~{output_table} \
                --barcode-umi-counts ~{barcode_umi_counts_file} \
                --knee-plot-pdf ~{knee_plot_pdf_file} \
                --cb-tag ~{cb_tag} \
                --umi-tag ~{umi_tag} \
                2>&1 | tee -a ~{log_file_name}
        else
            echo "Detected older cudll_bams_to_cb_umi.py CLI; generating summary outputs in WDL fallback" | tee -a ~{log_file_name}
            cudll_bams_to_cb_umi.py \
                --main-bam ~{main_bam} \
                ~{if defined(supp_bam) then "--supp-bam " + supp_bam else ""} \
                --output ~{output_table} \
                --cb-tag ~{cb_tag} \
                --umi-tag ~{umi_tag} \
                2>&1 | tee -a ~{log_file_name}

            python3 - "~{output_table}" "~{barcode_umi_counts_file}" "~{knee_plot_pdf_file}" "~{max_knee_plot_points}" <<'PY' 2>&1 | tee -a ~{log_file_name}
import collections
import gzip
import math
import sys


def open_text(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def pdf_escape(text):
    return text.replace("\\", "\\\\").replace("(", "\\(").replace(")", "\\)")


def select_plot_indices(num_points, max_points):
    if num_points <= max_points:
        return list(range(num_points))

    max_log_rank = math.log10(num_points) if num_points > 1 else 1.0
    selected = {0, num_points - 1}

    # Sample ranks across log space so the knee shape remains visible.
    for sample_idx in range(max_points):
        frac = sample_idx / (max_points - 1)
        rank = 10 ** (frac * max_log_rank)
        index = min(num_points - 1, max(0, int(round(rank)) - 1))
        selected.add(index)

    return sorted(selected)


def write_simple_pdf(pdf_path, barcode_counts, max_plot_points):
    width = 612
    height = 432
    left = 72
    right = width - 54
    bottom = 72
    top = height - 54
    plot_width = right - left
    plot_height = top - bottom

    cmds = [
        "0.0 0.0 0.0 RG",
        "1 w",
        f"{left} {bottom} m {left} {top} l S",
        f"{left} {bottom} m {right} {bottom} l S",
        "BT /F1 16 Tf 180 392 Td (Barcode UMI Knee Plot) Tj ET",
        "BT /F1 10 Tf 230 40 Td (Barcode rank) Tj ET",
        "q 1 0 0 1 20 220 cm 90 0 0 90 0 0 cm BT /F1 10 Tf 0 0 Td (Unique UMI count) Tj ET Q",
    ]

    if barcode_counts:
        counts_only = [count for _, count in barcode_counts]
        max_rank = len(counts_only)
        max_count = max(counts_only)
        min_positive = min(count for count in counts_only if count > 0)

        min_log_rank = 0.0
        max_log_rank = math.log10(max_rank) if max_rank > 1 else 1.0
        min_log_count = math.log10(min_positive)
        max_log_count = math.log10(max_count) if max_count > min_positive else min_log_count + 1.0

        plot_indices = select_plot_indices(len(counts_only), max_plot_points)
        points = []
        for index in plot_indices:
            idx = index + 1
            count = counts_only[index]
            x_frac = 0.0 if max_log_rank == min_log_rank else (math.log10(idx) - min_log_rank) / (max_log_rank - min_log_rank)
            y_frac = 0.5 if max_log_count == min_log_count else (math.log10(count) - min_log_count) / (max_log_count - min_log_count)
            x = left + x_frac * plot_width
            y = bottom + y_frac * plot_height
            points.append((x, y))

        if points:
            x0, y0 = points[0]
            cmds.append("0.1 0.3 0.7 RG")
            cmds.append("1.5 w")
            cmds.append(f"{x0:.2f} {y0:.2f} m")
            for x, y in points[1:]:
                cmds.append(f"{x:.2f} {y:.2f} l")
            cmds.append("S")

        top_count = counts_only[0]
        last_count = counts_only[-1]
        cmds.extend([
            "0.0 0.0 0.0 RG",
            f"BT /F1 9 Tf 82 {top + 6} Td ({pdf_escape(str(top_count))}) Tj ET",
            f"BT /F1 9 Tf {right - 24:.0f} 58 Td ({pdf_escape(str(max_rank))}) Tj ET",
            f"BT /F1 9 Tf 82 58 Td ({pdf_escape(str(last_count))}) Tj ET",
        ])
    else:
        cmds.append("BT /F1 12 Tf 180 220 Td (No valid barcode/UMI pairs found) Tj ET")

    stream = "\n".join(cmds) + "\n"
    objects = [
        "<< /Type /Catalog /Pages 2 0 R >>",
        "<< /Type /Pages /Kids [3 0 R] /Count 1 >>",
        f"<< /Type /Page /Parent 2 0 R /MediaBox [0 0 {width} {height}] /Resources << /Font << /F1 4 0 R >> >> /Contents 5 0 R >>",
        "<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica >>",
        f"<< /Length {len(stream.encode('latin-1'))} >>\nstream\n{stream}endstream",
    ]

    with open(pdf_path, "wb") as handle:
        handle.write(b"%PDF-1.4\n")
        offsets = [0]
        for idx, obj in enumerate(objects, start=1):
            offsets.append(handle.tell())
            handle.write(f"{idx} 0 obj\n".encode("latin-1"))
            handle.write(obj.encode("latin-1"))
            handle.write(b"\nendobj\n")
        xref_start = handle.tell()
        handle.write(f"xref\n0 {len(objects) + 1}\n".encode("latin-1"))
        handle.write(b"0000000000 65535 f \n")
        for offset in offsets[1:]:
            handle.write(f"{offset:010d} 00000 n \n".encode("latin-1"))
        handle.write(
            (
                f"trailer\n<< /Size {len(objects) + 1} /Root 1 0 R >>\n"
                f"startxref\n{xref_start}\n%%EOF\n"
            ).encode("latin-1")
        )


table_path, counts_path, pdf_path, max_plot_points_raw = sys.argv[1:]
max_plot_points = max(2, int(max_plot_points_raw))
unique_umis = collections.defaultdict(set)
with open_text(table_path) as handle:
    header = handle.readline().rstrip("\n")
    if header != "read_name\tcell_barcode\tUMI":
        raise RuntimeError(f"Unexpected CB/UMI header: {header}")
    for line in handle:
        line = line.rstrip("\n")
        if not line:
            continue
        _read_name, cell_barcode, umi = line.split("\t")
        if cell_barcode == "NA" or umi == "NA":
            continue
        unique_umis[cell_barcode].add(umi)

barcode_counts = sorted(
    ((cell_barcode, len(umis)) for cell_barcode, umis in unique_umis.items()),
    key=lambda item: (-item[1], item[0]),
)

with open(counts_path, "wt") as handle:
    handle.write("cell_barcode\tumi_count\n")
    for cell_barcode, umi_count in barcode_counts:
        handle.write(f"{cell_barcode}\t{umi_count}\n")

write_simple_pdf(pdf_path, barcode_counts, max_plot_points)
print(f"Wrote fallback barcode counts: {counts_path}")
print(f"Wrote fallback knee plot PDF: {pdf_path}")
PY
        fi

        echo "---" | tee -a ~{log_file_name}
        echo "Conversion complete" | tee -a ~{log_file_name}

        # Report output file size and line count
        ls -lh ~{output_table} | tee -a ~{log_file_name}
        ls -lh ~{barcode_umi_counts_file} | tee -a ~{log_file_name}
        ls -lh ~{knee_plot_pdf_file} | tee -a ~{log_file_name}
        echo "Line count:" | tee -a ~{log_file_name}
        zcat ~{output_table} | wc -l | tee -a ~{log_file_name}
    >>>

    output {
        File cb_umi_table = output_table
        File barcode_umi_counts = barcode_umi_counts_file
        File barcode_knee_plot_pdf = knee_plot_pdf_file
        File log_file = log_file_name
    }

    runtime {
        docker: docker
        cpu: cpu
        memory: "~{memory_gb} GB"
        disks: "local-disk ~{disk_gb} HDD"
        preemptible: 2
    }
}
