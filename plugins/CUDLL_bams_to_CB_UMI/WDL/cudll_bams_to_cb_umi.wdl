version 1.0

workflow CUDLLBamsToCbUmi {
    input {
        File cudll_main_bam
        File? cudll_supp_bam
        String sample_name
        String cb_tag = "CB"
        String umi_tag = "XM"
        Boolean generate_summary_outputs = true
        Int max_knee_plot_points = 5000
        Int sort_threads = 2
        String sort_memory_per_thread = "2G"
        String docker = "trinityctat/cudll-to-cb-umi"
        Int cpu = 8
        Int memory_gb = 32
        Int min_disk_gb = 100
        Int extra_disk_gb = 50
        Float disk_scale_factor = 6.0
        Int preemptible = 2
    }

    Float total_input_bam_gb = size(cudll_main_bam, "GB") + (if defined(cudll_supp_bam) then size(cudll_supp_bam, "GB") else 0.0)
    Int estimated_disk_gb = ceil(total_input_bam_gb * disk_scale_factor) + extra_disk_gb
    Int requested_disk_gb = if estimated_disk_gb > min_disk_gb then estimated_disk_gb else min_disk_gb

    if (generate_summary_outputs) {
        call BamsToCbUmiFull {
            input:
                main_bam = cudll_main_bam,
                supp_bam = cudll_supp_bam,
                sample_name = sample_name,
                cb_tag = cb_tag,
                umi_tag = umi_tag,
                max_knee_plot_points = max_knee_plot_points,
                sort_threads = sort_threads,
                sort_memory_per_thread = sort_memory_per_thread,
                docker = docker,
                cpu = cpu,
                memory_gb = memory_gb,
                disk_gb = requested_disk_gb,
                preemptible = preemptible
        }
    }

    if (!generate_summary_outputs) {
        call BamsToCbUmiPrimaryOnly {
            input:
                main_bam = cudll_main_bam,
                supp_bam = cudll_supp_bam,
                sample_name = sample_name,
                cb_tag = cb_tag,
                umi_tag = umi_tag,
                sort_threads = sort_threads,
                sort_memory_per_thread = sort_memory_per_thread,
                docker = docker,
                cpu = cpu,
                memory_gb = memory_gb,
                disk_gb = requested_disk_gb,
                preemptible = preemptible
        }
    }

    output {
        File cb_umi_table = select_first([BamsToCbUmiFull.cb_umi_table, BamsToCbUmiPrimaryOnly.cb_umi_table])
        File? barcode_umi_counts = BamsToCbUmiFull.barcode_umi_counts
        File? barcode_knee_plot_pdf = BamsToCbUmiFull.barcode_knee_plot_pdf
        File log_file = select_first([BamsToCbUmiFull.log_file, BamsToCbUmiPrimaryOnly.log_file])
    }
}

task BamsToCbUmiFull {
    input {
        File main_bam
        File? supp_bam
        String sample_name
        String cb_tag
        String umi_tag
        Int max_knee_plot_points
        Int sort_threads
        String sort_memory_per_thread
        String docker
        Int cpu
        Int memory_gb
        Int disk_gb
        Int preemptible
    }

    String output_table = "~{sample_name}.cb_umi.tsv.gz"
    String barcode_umi_counts_file = "~{sample_name}.barcode_umi_counts.tsv.gz"
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
        echo "Sort threads: ~{sort_threads}" | tee -a ~{log_file_name}
        echo "Sort memory per thread: ~{sort_memory_per_thread}" | tee -a ~{log_file_name}
        echo "Requested disk GB: ~{disk_gb}" | tee -a ~{log_file_name}
        echo "Preemptible attempts: ~{preemptible}" | tee -a ~{log_file_name}
        echo "---" | tee -a ~{log_file_name}

        cudll_bams_to_cb_umi.py \
            --main-bam ~{main_bam} \
            ~{if defined(supp_bam) then "--supp-bam " + supp_bam else ""} \
            --output ~{output_table} \
            --barcode-umi-counts ~{barcode_umi_counts_file} \
            --knee-plot-pdf ~{knee_plot_pdf_file} \
            --max-knee-plot-points ~{max_knee_plot_points} \
            --sort-threads ~{sort_threads} \
            --sort-memory-per-thread ~{sort_memory_per_thread} \
            --cb-tag ~{cb_tag} \
            --umi-tag ~{umi_tag} \
            2>&1 | tee -a ~{log_file_name}

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
        preemptible: preemptible
    }
}

task BamsToCbUmiPrimaryOnly {
    input {
        File main_bam
        File? supp_bam
        String sample_name
        String cb_tag
        String umi_tag
        Int sort_threads
        String sort_memory_per_thread
        String docker
        Int cpu
        Int memory_gb
        Int disk_gb
        Int preemptible
    }

    String output_table = "~{sample_name}.cb_umi.tsv.gz"
    String log_file_name = "~{sample_name}.cudll_cb_umi.log"

    command <<<
        set -e
        set -o pipefail

        echo "Starting CUDLL BAM to CB/UMI table conversion" | tee ~{log_file_name}
        echo "Sample: ~{sample_name}" | tee -a ~{log_file_name}
        echo "Main BAM: ~{main_bam}" | tee -a ~{log_file_name}
        echo "Supp BAM: ~{if defined(supp_bam) then supp_bam else "none"}" | tee -a ~{log_file_name}
        echo "Output: ~{output_table}" | tee -a ~{log_file_name}
        echo "Barcode UMI counts: disabled" | tee -a ~{log_file_name}
        echo "Knee plot PDF: disabled" | tee -a ~{log_file_name}
        echo "CB Tag: ~{cb_tag}" | tee -a ~{log_file_name}
        echo "UMI Tag: ~{umi_tag}" | tee -a ~{log_file_name}
        echo "Sort threads: ~{sort_threads}" | tee -a ~{log_file_name}
        echo "Sort memory per thread: ~{sort_memory_per_thread}" | tee -a ~{log_file_name}
        echo "Requested disk GB: ~{disk_gb}" | tee -a ~{log_file_name}
        echo "Preemptible attempts: ~{preemptible}" | tee -a ~{log_file_name}
        echo "---" | tee -a ~{log_file_name}

        cudll_bams_to_cb_umi.py \
            --main-bam ~{main_bam} \
            ~{if defined(supp_bam) then "--supp-bam " + supp_bam else ""} \
            --output ~{output_table} \
            --sort-threads ~{sort_threads} \
            --sort-memory-per-thread ~{sort_memory_per_thread} \
            --cb-tag ~{cb_tag} \
            --umi-tag ~{umi_tag} \
            2>&1 | tee -a ~{log_file_name}

        echo "---" | tee -a ~{log_file_name}
        echo "Conversion complete" | tee -a ~{log_file_name}

        ls -lh ~{output_table} | tee -a ~{log_file_name}
        echo "Line count:" | tee -a ~{log_file_name}
        zcat ~{output_table} | wc -l | tee -a ~{log_file_name}
    >>>

    output {
        File cb_umi_table = output_table
        File log_file = log_file_name
    }

    runtime {
        docker: docker
        cpu: cpu
        memory: "~{memory_gb} GB"
        disks: "local-disk ~{disk_gb} HDD"
        preemptible: preemptible
    }
}
