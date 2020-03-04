import "https://api.firecloud.org/ga4gh/v1/tools/CTAT:star_fusion_tasks/versions/1/plain-WDL/descriptor" as star_fusion_tasks

workflow fusion_inspector_workflow {
    File fusion_predictions
    File left_fq
    File? right_fq
    String sample_id
    # either a gs:// URL to a tar.gz2 file or an acronym
    String genome

    # cpus and memory defaults are read from acronym file by default.
    Int? num_cpu
    String? memory

    String? docker = "trinityctat/starfusion:1.8.1b"

    Float? extra_disk_space = 10
    Float? fastq_disk_space_multiplier = 3.25
    Float? genome_disk_space_multiplier = 2.5

    String? additional_flags

    String? acronym_file = "gs://regev-lab/resources/ctat/star_fusion/index.json"

    Int? preemptible = 2
    Boolean? use_ssd = true
    String config_docker = "continuumio/miniconda3:4.6.14"
    String? zones = "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"

    call star_fusion_tasks.star_fusion_config {
        input:
            genome=genome,
            acronym_file=acronym_file,
            docker=config_docker,
            cpus=num_cpu,
            memory=memory,
            preemptible = preemptible
    }

    call fusion_inspector {
        input:
            fusion_predictions=fusion_predictions,
            left_fq=left_fq,
            right_fq=right_fq,
            sample_id=sample_id,
            additional_flags=additional_flags,
            extra_disk_space=extra_disk_space,
            fastq_disk_space_multiplier=fastq_disk_space_multiplier,
            genome_disk_space_multiplier=genome_disk_space_multiplier,
            zones = zones,
            preemptible = preemptible,
            docker = docker,
            genome=star_fusion_config.star_genome,
            cpu=star_fusion_config.star_cpus_output,
            memory=star_fusion_config.star_memory_output,
            use_ssd=use_ssd
     }


    output {
        File fusion_inspector_inspect_web  = fusion_inspector.fusion_inspector_inspect_web
        File fusion_inspector_inspect_fusions_abridged = fusion_inspector.fusion_inspector_inspect_fusions_abridged
    }
}


task fusion_inspector {
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
    String? fusion_inspector
    String? additional_flags
    Boolean use_ssd

    command {
        set -e

        mkdir -p ${sample_id}
        mkdir -p genome_dir

        pbzip2 -dc ${genome} | tar x -C genome_dir --strip-components 1

       /usr/local/src/STAR-Fusion/FusionInspector/FusionInspector \
        --fusions ${fusion_predictions} \
        --genome_lib_dir `pwd`/genome_dir/ctat_genome_lib_build_dir \
        -O ${sample_id} \
        --CPU ${cpu} \
        --left_fq ${left_fq} \
        ${"--right_fq " + right_fq} \
        --out_prefix finspector \
        --vis \
        ${"" + additional_flags}
    }

    output {
        File fusion_inspector_inspect_web  = "${sample_id}/finspector.fusion_inspector_web.html"
        File fusion_inspector_inspect_fusions_abridged = "${sample_id}/finspector.FusionInspector.fusions.abridged.tsv"
    }

    runtime {
        docker: "${docker}"
        zones: zones
        disks: "local-disk " + ceil((fastq_disk_space_multiplier * (size(left_fq,  "GB") + size(right_fq,  "GB"))) +size(genome, "GB")*genome_disk_space_multiplier + extra_disk_space)+ " " + (if use_ssd then "SSD" else "HDD")
        memory :"${memory}"
        preemptible: "${preemptible}"
        cpu:"${cpu}"
    }
}
