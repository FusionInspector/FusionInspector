version 1.0

workflow CUDLLBamsToCbUmi {
    input {
        File cudll_main_bam
        File? cudll_supp_bam
        String sample_name
        String cb_tag = "CB"
        String umi_tag = "XM"
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
        echo "---" | tee -a ~{log_file_name}

        # Run the conversion script
        cudll_bams_to_cb_umi.py \
            --main-bam ~{main_bam} \
            ~{if defined(supp_bam) then "--supp-bam " + supp_bam else ""} \
            --output ~{output_table} \
            --barcode-umi-counts ~{barcode_umi_counts_file} \
            --knee-plot-pdf ~{knee_plot_pdf_file} \
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
        preemptible: 2
    }
}
