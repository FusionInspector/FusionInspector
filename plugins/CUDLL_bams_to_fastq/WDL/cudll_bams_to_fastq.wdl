version 1.0

workflow CUDLLBamsToFastq {
    input {
        File cudll_main_bam
        File cudll_supp_bam
        String sample_name
        String docker = "trinityctat/cudll-to-fastq"
        Int cpu = 4
        Int memory_gb = 16
        Int disk_gb = 100
    }

    call BamsToFastq {
        input:
            main_bam = cudll_main_bam,
            supp_bam = cudll_supp_bam,
            sample_name = sample_name,
            docker = docker,
            cpu = cpu,
            memory_gb = memory_gb,
            disk_gb = disk_gb
    }

    output {
        File merged_fastq = BamsToFastq.merged_fastq
        File log_file = BamsToFastq.log_file
    }
}

task BamsToFastq {
    input {
        File main_bam
        File supp_bam
        String sample_name
        String docker
        Int cpu
        Int memory_gb
        Int disk_gb
    }

    String output_fastq = "~{sample_name}.merged.fastq.gz"
    String log_file_name = "~{sample_name}.cudll_conversion.log"

    command <<<
        set -e
        set -o pipefail

        echo "Starting CUDLL BAM to FASTQ conversion" | tee ~{log_file_name}
        echo "Sample: ~{sample_name}" | tee -a ~{log_file_name}
        echo "Main BAM: ~{main_bam}" | tee -a ~{log_file_name}
        echo "Supp BAM: ~{supp_bam}" | tee -a ~{log_file_name}
        echo "Output: ~{output_fastq}" | tee -a ~{log_file_name}
        echo "---" | tee -a ~{log_file_name}

        # Run the conversion script
        cudll_bams_to_fastq.py \
            --main-bam ~{main_bam} \
            --supp-bam ~{supp_bam} \
            --output ~{output_fastq} \
            2>&1 | tee -a ~{log_file_name}

        echo "---" | tee -a ~{log_file_name}
        echo "Conversion complete" | tee -a ~{log_file_name}

        # Report output file size
        ls -lh ~{output_fastq} | tee -a ~{log_file_name}
    >>>

    output {
        File merged_fastq = output_fastq
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
