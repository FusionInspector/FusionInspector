version 1.0


workflow fusion_inspector_workflow {

  input {
    String? docker = "trinityctat/fusioninspector:latest"
    Int? num_cpu
    String? memory
    Boolean? use_ssd = true
    String genome
    Float? genome_disk_space_multiplier = 2.5
    Float? fastq_disk_space_multiplier = 3.25
    File fusion_predictions
    
    String? additional_flags
    Int? preemptible = 2
    File? right_fq
    Float? extra_disk_space = 10
    String? zones = "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
    String sample_id


    File? rnaseq_aligned_bam
    File? fastq_pair_tar_gz
    File? left_fq
    File? right_fq
        
  }

  
   if (defined(rnaseq_aligned_bam)) {
      
      call CTAT_BAM_TO_FASTQ {
            input:
              input_bam=rnaseq_aligned_bam,
              sample_name=sample_name,
              Docker=Docker
        }
    }

    if (defined(fastq_pair_tar_gz)) {

      call CTAT_UNTAR_FASTQS {
        input:
          fastq_pair_tar_gz=fastq_pair_tar_gz,
          Docker=Docker
      }
    }

    
    File? left_fq_use = select_first([left_fq, CTAT_UNTAR_FASTQS.left_fq, CTAT_BAM_TO_FASTQ.left_fq])
    File? right_fq_use = select_first([right_fq, CTAT_UNTAR_FASTQS.right_fq, CTAT_BAM_TO_FASTQ.right_fq])

  
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

task CTAT_BAM_TO_FASTQ {

    input {
      File input_bam
      String sample_name
      String Docker
    }
    
    command {

    set -e

    # initial potential cleanup of read names in the bam file
    /usr/local/bin/sam_readname_cleaner.py ${input_bam} ${input_bam}.cleaned.bam


    # revert aligned bam
    java -Xmx1000m -jar /usr/local/src/picard.jar \
        RevertSam \
        INPUT=${input_bam}.cleaned.bam \
        OUTPUT_BY_READGROUP=false \
        VALIDATION_STRINGENCY=SILENT \
        SORT_ORDER=queryname \
        OUTPUT=${sample_name}.reverted.bam 


    # bam to fastq
    java -jar /usr/local/src/picard.jar \
        SamToFastq I=${sample_name}.reverted.bam \
        F=${sample_name}_1.fastq F2=${sample_name}_2.fastq \
        INTERLEAVE=false NON_PF=true \
        CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2

   }
    
    output {
      File left_fq="${sample_name}_1.fastq"
      File right_fq="${sample_name}_2.fastq"
    }
    
    runtime {
            docker: Docker
            disks: "local-disk 500 SSD"
            memory: "20G"
            cpu: "16"
            preemptible: 0
            maxRetries: 3
    }
    
  }


  
task CTAT_UNTAR_FASTQS {

  input {
    File fastq_pair_tar_gz
    String Docker
  }

  command {

     set -e
    
     # untar the fq pair
     tar xvf ${fastq_pair_tar_gz}
  }

  output {
    File left_fq = select_first(glob("*_1.fastq"))
    File right_fq = select_first(glob("*_2.fastq"))
  }

  runtime {
            docker: Docker
            disks: "local-disk 500 SSD"
            memory: "10G"
            cpu: "4"
            preemptible: 0
            maxRetries: 3
    }
  }


  
