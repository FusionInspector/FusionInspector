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
        Int sort_threads = 1
        String sort_memory_per_thread = "1G"
        String docker = "trinityctat/cudll-to-cb-umi"
        Int cpu = 4
        Int memory_gb = 16
        Int disk_gb = 100
        Int preemptible = 2
    }

    call BamsToCbUmi {
        input:
            main_bam = cudll_main_bam,
            supp_bam = cudll_supp_bam,
            sample_name = sample_name,
            cb_tag = cb_tag,
            umi_tag = umi_tag,
            generate_summary_outputs = generate_summary_outputs,
            max_knee_plot_points = max_knee_plot_points,
            sort_threads = sort_threads,
            sort_memory_per_thread = sort_memory_per_thread,
            docker = docker,
            cpu = cpu,
            memory_gb = memory_gb,
            disk_gb = disk_gb,
            preemptible = preemptible
    }

    output {
        File cb_umi_table = BamsToCbUmi.cb_umi_table
        File? barcode_umi_counts = BamsToCbUmi.barcode_umi_counts
        File? barcode_knee_plot_pdf = BamsToCbUmi.barcode_knee_plot_pdf
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
        Boolean generate_summary_outputs
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
        echo "Generate summary outputs: ~{generate_summary_outputs}" | tee -a ~{log_file_name}
        echo "Max knee plot points: ~{max_knee_plot_points}" | tee -a ~{log_file_name}
        echo "Sort threads: ~{sort_threads}" | tee -a ~{log_file_name}
        echo "Sort memory per thread: ~{sort_memory_per_thread}" | tee -a ~{log_file_name}
        echo "Preemptible attempts: ~{preemptible}" | tee -a ~{log_file_name}
        echo "---" | tee -a ~{log_file_name}

        cudll_bams_to_cb_umi.py \
            --main-bam ~{main_bam} \
            ~{if defined(supp_bam) then "--supp-bam " + supp_bam else ""} \
            --output ~{output_table} \
            ~{if generate_summary_outputs then "--barcode-umi-counts " + barcode_umi_counts_file else ""} \
            ~{if generate_summary_outputs then "--knee-plot-pdf " + knee_plot_pdf_file else ""} \
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
        if [ "~{generate_summary_outputs}" = "true" ]; then
            ls -lh ~{barcode_umi_counts_file} | tee -a ~{log_file_name}
            ls -lh ~{knee_plot_pdf_file} | tee -a ~{log_file_name}
        fi
        echo "Line count:" | tee -a ~{log_file_name}
        zcat ~{output_table} | wc -l | tee -a ~{log_file_name}
    >>>

    output {
        File cb_umi_table = output_table
        File? barcode_umi_counts = if generate_summary_outputs then barcode_umi_counts_file else None
        File? barcode_knee_plot_pdf = if generate_summary_outputs then knee_plot_pdf_file else None
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
