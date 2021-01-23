version 1.0


workflow fusion_inspector_workflow {

  input {
    String? docker = "trinityctat/fusioninspector:latest"
    Int? num_cpu = 10
    String? memory = "50G"
    Boolean? use_ssd = true
    String genome
    Float? genome_disk_space_multiplier = 2.5
    Float? fastq_disk_space_multiplier = 3.25
    File fusion_predictions
    String sample_name
    
    String? additional_flags
    Int? preemptible = 2
    Int? maxRetries = 2
    File? right_fq
    Float? extra_disk_space = 10
    String? zones = "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
   


    File? rnaseq_aligned_bam
    File? fastq_pair_tar_gz
    File? left_fq
    File? right_fq
        
  }

  
   if (defined(rnaseq_aligned_bam)) {
      
      call CTAT_BAM_TO_FASTQ {
            input:
              input_bam=select_first([rnaseq_aligned_bam]),
              sample_name=sample_name,
              Docker=select_first([docker]),
              maxRetries=select_first([maxRetries]),
              preemptible=select_first([preemptible])
        }
    }

    if (defined(fastq_pair_tar_gz)) {

      call CTAT_UNTAR_FASTQS {
        input:
          fastq_pair_tar_gz=select_first([fastq_pair_tar_gz]),
          Docker=select_first([docker]),
          maxRetries=select_first([maxRetries]),
          preemptible=select_first([preemptible])
          
      }
    }

    
    File left_fq_use = select_first([left_fq, CTAT_UNTAR_FASTQS.left_fq, CTAT_BAM_TO_FASTQ.left_fq])
    File right_fq_use = select_first([right_fq, CTAT_UNTAR_FASTQS.right_fq, CTAT_BAM_TO_FASTQ.right_fq])

  
  call fusion_inspector {
    input:
      fusion_predictions = fusion_predictions,
      genome = genome,
      sample_id = sample_name,
      left_fq = left_fq_use,
      right_fq = right_fq_use,
      zones = select_first([zones]),
      preemptible = select_first([preemptible]),
      docker = select_first([docker]),
      cpu =  select_first([num_cpu]),
      memory = select_first([memory]),
      extra_disk_space = select_first([extra_disk_space]),
      fastq_disk_space_multiplier = select_first([fastq_disk_space_multiplier]),
      genome_disk_space_multiplier = select_first([genome_disk_space_multiplier]),
      additional_flags = additional_flags,
      use_ssd = select_first([use_ssd]),
      maxRetries = select_first([maxRetries])
  }


  output {
    File fusion_inspector_web = fusion_inspector.fusion_inspector_web
    File fusion_inspector_fusions_abridged = fusion_inspector.fusion_inspector_fusions_abridged
    
    File fusion_inspector_contigs_fa = fusion_inspector.fusion_inspector_contigs_fa 
    File fusion_inspector_consolidated_bam = fusion_inspector.fusion_inspector_consolidated_bam
    File fusion_inspector_consolidated_bam_bai = fusion_inspector.fusion_inspector_consolidated_bam_bai
    File fusion_inspector_junction_reads_bam = fusion_inspector.fusion_inspector_junction_reads_bam
    File fusion_inspector_junction_reads_bam_bai = fusion_inspector.fusion_inspector_junction_reads_bam_bai 
    File fusion_inspector_spanning_reads_bam = fusion_inspector.fusion_inspector_spanning_reads_bam
    File fusion_inspector_spanning_reads_bam_bai= fusion_inspector.fusion_inspector_spanning_reads_bam_bai
    
    File fusion_inspector_genes_bed = fusion_inspector.fusion_inspector_genes_bed
    File fusion_inspector_pfam_bed = fusion_inspector.fusion_inspector_pfam_bed
    File fusion_inspector_seqsimilar_bed = fusion_inspector.fusion_inspector_seqsimilar_bed
    
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
    Int maxRetries
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
    File fusion_inspector_web = "${sample_id}/${sample_id}.fusion_inspector_web.html"
    File fusion_inspector_fusions_abridged = "${sample_id}/${sample_id}.FusionInspector.fusions.abridged.tsv"
    
    File fusion_inspector_contigs_fa = "${sample_id}/${sample_id}.fa"
    File fusion_inspector_consolidated_bam = "${sample_id}/${sample_id}.consolidated.bam"
    File fusion_inspector_consolidated_bam_bai = "${sample_id}/${sample_id}.consolidated.bam.bai"
    File fusion_inspector_junction_reads_bam = "${sample_id}/${sample_id}.junction_reads.bam"
    File fusion_inspector_junction_reads_bam_bai = "${sample_id}/${sample_id}.junction_reads.bam.bai"
    File fusion_inspector_spanning_reads_bam = "${sample_id}/${sample_id}.spanning_reads.bam"
    File fusion_inspector_spanning_reads_bam_bai = "${sample_id}/${sample_id}.spanning_reads.bam.bai"
    
    File fusion_inspector_genes_bed = "${sample_id}/${sample_id}.bed"
    File fusion_inspector_pfam_bed = "${sample_id}/${sample_id}.igv.Pfam.bed"
    File fusion_inspector_seqsimilar_bed = "${sample_id}/${sample_id}.igv.seqsimilar.bed"
    
  }
  command <<<

        set -e

        mkdir -p ~{sample_id}
        mkdir -p genome_dir

        pbzip2 -dc ~{genome} | tar x -C genome_dir --strip-components 1

        FusionInspector \
        --fusions ~{fusion_predictions} \
        --genome_lib_dir `pwd`/genome_dir \
        -O ~{sample_id} \
        --CPU ~{cpu} \
        --left_fq ~{left_fq} \
        ~{"--right_fq " + right_fq} \
        --out_prefix ~{sample_id} \
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
    maxRetries: "${maxRetries}"
  }

}

task CTAT_BAM_TO_FASTQ {

    input {
      File? input_bam
      String sample_name
      String Docker
      Int preemptible
      Int maxRetries
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
            preemptible: "${preemptible}"
            maxRetries: "${maxRetries}"
    }
    
}


  
task CTAT_UNTAR_FASTQS {

  input {
    File fastq_pair_tar_gz
    String Docker
    Int preemptible
    Int maxRetries
  }

  command {

     set -e
    
     # untar the fq pair
     tar xvf ${fastq_pair_tar_gz}
  }

  output {
    File left_fq = glob("*_1.fastq*")[0]
    File right_fq = glob("*_2.fastq*")[0]
  }

  runtime {
            docker: Docker
            disks: "local-disk 500 SSD"
            memory: "10G"
            cpu: "4"
            preemptible: "${preemptible}"
            maxRetries: "${maxRetries}"
    }
  }


  
