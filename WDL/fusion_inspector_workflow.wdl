version 1.0
import "https://api.firecloud.org/ga4gh/v1/tools/CTAT:star_fusion_tasks/versions/2/plain-WDL/descriptor" as star_fusion_tasks

#import "./star_fusion_tasks.wdl" as star_fusion_tasks

workflow fusion_inspector_workflow {
  input {
    String? docker = "trinityctat/starfusion:1.8.1b"
    Int? num_cpu
    String config_docker = "continuumio/miniconda3:4.6.14"
    String? memory
    Boolean? use_ssd = true
    String genome
    Float? genome_disk_space_multiplier = 2.5
    String? acronym_file = "gs://regev-lab/resources/ctat/star_fusion/index.json"
    Float? fastq_disk_space_multiplier = 3.25
    File fusion_predictions
    File left_fq
    String? additional_flags
    Int? preemptible = 2
    File? right_fq
    Float? extra_disk_space = 10
    String? zones = "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
    String sample_id
  }

   call star_fusion_tasks.star_fusion_config as star_fusion_config {
      input:
        genome = genome,
        acronym_file = acronym_file,
        cpus = num_cpu,
        memory = memory,
        docker = config_docker,
        preemptible = preemptible
    }

  call fusion_inspector {
    input:
      fusion_predictions = fusion_predictions,
      genome = star_fusion_config.star_genome,
      sample_id = sample_id,
      left_fq = left_fq,
      right_fq = right_fq,
      zones = zones,
      preemptible = preemptible,
      docker = docker,
      cpu = star_fusion_config.star_cpus_output,
      memory = star_fusion_config.star_memory_output,
      extra_disk_space = extra_disk_space,
      fastq_disk_space_multiplier = fastq_disk_space_multiplier,
      genome_disk_space_multiplier = genome_disk_space_multiplier,
      additional_flags = additional_flags,
      use_ssd = use_ssd
  }


  output {
    File fusion_inspector_inspect_web = fusion_inspector.fusion_inspector_inspect_web
    File fusion_inspector_inspect_fusions_abridged = fusion_inspector.fusion_inspector_inspect_fusions_abridged
  }
}
task fusion_inspector {
  input {
    File fusion_predictions
    File genome
    String sample_id
    File left_fq
    File right_fq
    String zones
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


  output {
    File fusion_inspector_inspect_web = "${sample_id}/finspector.fusion_inspector_web.html"
    File fusion_inspector_inspect_fusions_abridged = "${sample_id}/finspector.FusionInspector.fusions.abridged.tsv"
  }
  command <<<

        set -e

        mkdir -p ~{sample_id}
        mkdir -p genome_dir

        pbzip2 -dc ~{genome} | tar x -C genome_dir --strip-components 1

       /usr/local/src/STAR-Fusion/FusionInspector/FusionInspector \
        --fusions ~{fusion_predictions} \
        --genome_lib_dir `pwd`/genome_dir/ctat_genome_lib_build_dir \
        -O ~{sample_id} \
        --CPU ~{cpu} \
        --left_fq ~{left_fq} \
        ~{"--right_fq " + right_fq} \
        --out_prefix finspector \
        --vis \
        ~{"" + additional_flags}

  >>>
  runtime {
    preemptible: "${preemptible}"
    disks: "local-disk " + ceil((fastq_disk_space_multiplier * (size(left_fq, "GB") + size(right_fq, "GB"))) + size(genome, "GB") * genome_disk_space_multiplier + extra_disk_space) + " " + (if use_ssd then "SSD" else "HDD")
    docker: "${docker}"
    cpu: "${cpu}"
    zones: zones
    memory: "${memory}"
  }

}

