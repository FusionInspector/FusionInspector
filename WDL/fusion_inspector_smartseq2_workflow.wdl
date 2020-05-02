version 1.0
import "https://api.firecloud.org/ga4gh/v1/tools/jgould:star_fusion_tasks/versions/4/plain-WDL/descriptor" as star_fusion_tasks

#import "/Users/jgould/git/STAR-Fusion/WDL/star_fusion_tasks.wdl" as star_fusion_tasks

workflow fusion_inspector_smartseq2_workflow {
  input {
    String? docker
    String? config_docker
    String? gatk_docker
    Boolean? use_ssd
    Int? preemptible
    Int? num_cpu
    String? acronym_file
    String? memory
    String genome
    Float genome_disk_space_multiplier = 2.5
    Float fastq_disk_space_multiplier = 3.25
    Float merge_bam_disk_space_multiplier = 3.25
    Float merge_bam_extra_disk_space = 2
    File fusion_predictions

    String? additional_flags
    File sample_sheet # csv with Cell,Read1,Read2
    Int? cells_per_job
    Float extra_disk_space = 10
    String sample_id
  }

  Int preemptible_or_default = select_first([preemptible, 2])
  Int cells_per_job_or_default = select_first([cells_per_job, 24])
  Boolean use_ssd_or_default = select_first([use_ssd, true])

  String docker_or_default = select_first([docker, "bhaastestdockers/fusioninspector:2.4.0-dev"])
  String acronym_file_or_default = select_first([acronym_file, "gs://regev-lab/resources/ctat/star_fusion/index.json"])
  String config_docker_or_default = select_first([config_docker, "continuumio/anaconda3:2020.02"])
  String gatk_docker_or_default = select_first([gatk_docker, 'us.gcr.io/broad-gatk/gatk:4.1.7.0'])

  call star_fusion_tasks.star_fusion_config as star_fusion_config {
      input:
        genome = genome,
        acronym_file = acronym_file_or_default,
        cpus = num_cpu,
        memory = memory,
        docker = config_docker_or_default,
        preemptible = preemptible_or_default
  }

   call star_fusion_tasks.split_sample_sheet {
        input:
            sample_sheet=sample_sheet,
            cells_per_job=cells_per_job_or_default,
            docker=config_docker_or_default,
            preemptible = preemptible_or_default
   }

    scatter(idx in range(length(split_sample_sheet.split_output["read1"]))) {
        call fusion_inspector {
            input:
                left_fq=split_sample_sheet.split_output.read1[idx],
                right_fq=split_sample_sheet.split_output.read2[idx],
                cell_name=split_sample_sheet.split_output.cell[idx],
                fusion_predictions = fusion_predictions,
                genome = star_fusion_config.star_genome,
                sample_id = sample_id,
                preemptible = preemptible_or_default,
                docker = docker_or_default,
                cpu = star_fusion_config.star_cpus_output,
                memory = star_fusion_config.star_memory_output,
                extra_disk_space = extra_disk_space,
                fastq_disk_space_multiplier = fastq_disk_space_multiplier,
                genome_disk_space_multiplier = genome_disk_space_multiplier,
                additional_flags = additional_flags,
                use_ssd = use_ssd_or_default
        }
    }

    call star_fusion_tasks.concatentate as concatentate_fusion_inspector_inspect_fusions_samples_deconvolved {
        input:
             input_files=fusion_inspector.samples_deconvolved,
             output_name="fusions.samples_deconvolved.tsv",
             docker=config_docker_or_default,
             preemptible = preemptible_or_default
    }

    call merge_bam_files as merge_junction_reads {
        input:
            input_bams=fusion_inspector.junction_reads,
            output_name="finspector.junction_reads",
            preemptible=preemptible_or_default,
            docker=gatk_docker_or_default,
            memory="2G",
            extra_disk_space=merge_bam_extra_disk_space,
            disk_space_multiplier=merge_bam_disk_space_multiplier
    }

     call merge_bam_files as merge_spanning_reads {
        input:
            input_bams=fusion_inspector.spanning_reads,
            output_name="finspector.spanning_reads",
            preemptible=preemptible_or_default,
            docker=gatk_docker_or_default,
            memory="2G",
            extra_disk_space=merge_bam_extra_disk_space,
            disk_space_multiplier=merge_bam_disk_space_multiplier
        }


  output {
    File fa = fusion_inspector.fa[0]
    File fai = fusion_inspector.fai[0]
    File bed = fusion_inspector.bed[0]
    File gtf = fusion_inspector.gtf[0]
    File cytoband = fusion_inspector.cytoband[0]

    File fusions_samples_deconvolved = concatentate_fusion_inspector_inspect_fusions_samples_deconvolved.output_file
    File junction_reads = merge_junction_reads.bam
    File junction_reads_index = merge_junction_reads.bai
    File spanning_reads = merge_spanning_reads.bam
    File spanning_reads_index = merge_spanning_reads.bai
  }
}


task fusion_inspector {
  input {
    Array[File] left_fq
    Array[File] right_fq
    Array[String] cell_name
    File fusion_predictions
    File genome
    String sample_id
    Int preemptible
    String docker
    Int cpu
    String memory
    Float extra_disk_space
    Float fastq_disk_space_multiplier
    Float genome_disk_space_multiplier
    String? additional_flags
    Boolean use_ssd
  }

  command <<<

    set -e

    mkdir -p ~{sample_id}
    mkdir -p genome_dir

    python <<CODE
    left_fq = "~{sep=',' left_fq}".split(',')
    right_fq = "~{sep=',' right_fq}".split(',')
    cell_name = "~{sep=',' cell_name}".split(',')
    with open('samples_file.txt', 'wt') as f:
        for i in range(len(cell_name)):
            f.write(cell_name[i] + '\t')
            f.write(left_fq[i])
            if i < len(right_fq):
                f.write('\t')
                f.write(right_fq[i])
            f.write('\n')
    CODE

    pbzip2 -dc ~{genome} | tar x -C genome_dir --strip-components 1

    /usr/local/bin/FusionInspector \
    --fusions ~{fusion_predictions} \
    --genome_lib_dir `pwd`/genome_dir/ctat_genome_lib_build_dir \
    -O ~{sample_id} \
    --CPU ~{cpu} \
    --samples_file samples_file.txt \
    --out_prefix finspector \
    --vis \
    ~{"" + additional_flags}

  >>>

    output {
        File samples_deconvolved = "~{sample_id}/finspector.FusionInspector.fusions.samples_deconvolved.tsv"
        File junction_reads = "~{sample_id}/finspector.junction_reads.bam"
        File spanning_reads = "~{sample_id}/finspector.spanning_reads.bam"
        File fa = "~{sample_id}/finspector.fa"
        File fai = "~{sample_id}/finspector.fa.fai"
        File bed = "~{sample_id}/finspector.bed"
        File gtf = "~{sample_id}/finspector.gtf"
        File cytoband = "~{sample_id}/cytoBand.txt"


#          File? igv_pfam_bed = "~{sample_id}/finspector.igv.Pfam.bed"
#      File? fusion_junc_span = "~{sample_id}/finspector.igv.FusionJuncSpan"
#      File? junction_reads_bed = "~{sample_id}/finspector.junction_reads.bam.bed"
#      File? spanning_reads_bed = "~{sample_id}/finspector.spanning_reads.bam.bed"

   }

  runtime {
    preemptible: "~{preemptible}"
    disks: "local-disk " + ceil((fastq_disk_space_multiplier * (size(left_fq, "GB") + size(right_fq, "GB"))) + size(genome, "GB") * genome_disk_space_multiplier + extra_disk_space) + " " + (if use_ssd then "SSD" else "HDD")
    docker: "~{docker}"
    cpu: "~{cpu}"
    memory: "~{memory}"
  }
}


task merge_bam_files {
    input {
        Array[File] input_bams
        String output_name
        Int preemptible
        String docker
        String memory
        Float extra_disk_space
        Float disk_space_multiplier
    }

    command <<<
        mem=$(cat /proc/meminfo | grep MemAvailable | awk 'BEGIN { FS=" " } ; { print $2 }')
        mem="$(($mem/1000))"
        gatk --java-options "-Xmx$(echo $mem)m" MergeSamFiles \
              -CREATE_INDEX \
              -I=~{sep=' -I=' input_bams} \
              -O=~{output_name}.bam \
              -VALIDATION_STRINGENCY LENIENT
    >>>

    output {
        File bam = "~{output_name}.bam"
        File bai = "~{output_name}.bai"
    }

    runtime {
        preemptible: "~{preemptible}"
        disks: "local-disk " + ceil((disk_space_multiplier * (size(input_bams, "GB")))+ extra_disk_space) + " HDD"
        docker: "~{docker}"
        cpu: 1
        memory: "~{memory}"
    }
}
