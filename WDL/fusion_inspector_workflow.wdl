version 1.0


workflow fusion_inspector_workflow {

  input {

    String sample_id
    
    File genome_plug_n_play_tar_gz = "gs://mdl-ctat-genome-libs/__genome_libs_StarFv1.10/GRCh38_gencode_v22_CTAT_lib_Mar012021.plug-n-play.tar.gz"
    File left_fq
    File? right_fq
    File target_fusions_list    

    String docker = "trinityctat/fusioninspector:latest"
    
    String? additional_flags

    Int num_cpu = 12
    String memory = "50G"
    Boolean use_ssd = true
    Float genome_disk_space_multiplier = 2.5
    Float fastq_disk_space_multiplier = 3.25
    Int preemptible = 1
    Float extra_disk_space = 10    
        
  }
  
  call fusion_inspector {
    input:
      target_fusions_list = target_fusions_list,
      genome_plug_n_play_tar_gz = genome_plug_n_play_tar_gz,
      sample_id = sample_id,
      left_fq = left_fq,
      right_fq = right_fq,
    
      preemptible = preemptible,
      docker = docker,
      cpu = num_cpu,
      memory = memory,
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
    File target_fusions_list
    File genome_plug_n_play_tar_gz
    String sample_id
    File left_fq
    File? right_fq
   
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
    File fusion_inspector_inspect_web = "~{sample_id}.fusion_inspector_web.html"
    File fusion_inspector_inspect_fusions_abridged = "~{sample_id}.FusionInspector.fusions.abridged.tsv"
    File fusion_inspector_inspect_fusions = "~{sample_id}.FusionInspector.fusions.tsv"
    File fusion_inspector_IGV_inputs = "~{sample_id}.IGV_inputs.tar.gz"
  }

  command <<<

        set -ex

        mkdir -p ~{sample_id}
        mkdir -p genome_dir
     
        tar xf ~{genome_plug_n_play_tar_gz} -C genome_dir --strip-components 1
        
       FusionInspector \
        --fusions ~{target_fusions_list} \
        --genome_lib_dir `pwd`/genome_dir/ctat_genome_lib_build_dir \
        -O ~{sample_id} \
        --CPU ~{cpu} \
        --left_fq ~{left_fq} \
        ~{"--right_fq " + right_fq} \
        --vis \
        ~{"" + additional_flags}

        mv ~{sample_id}/IGV_inputs ~{sample_id}.IGV_inputs
        tar -zcvf ~{sample_id}.IGV_inputs.tar.gz ~{sample_id}.IGV_inputs

        mv ~{sample_id}/finspector.fusion_inspector_web.html  ~{sample_id}.fusion_inspector_web.html
        mv ~{sample_id}/finspector.FusionInspector.fusions.abridged.tsv ~{sample_id}.FusionInspector.fusions.abridged.tsv
        mv ~{sample_id}/finspector.FusionInspector.fusions.tsv ~{sample_id}.FusionInspector.fusions.abridged.tsv
        
    
  >>>

  
  runtime {
    preemptible: "${preemptible}"
    disks: "local-disk " + ceil((fastq_disk_space_multiplier * (size(left_fq, "GB") + size(right_fq, "GB"))) + size(genome_plug_n_play_tar_gz, "GB") * genome_disk_space_multiplier + extra_disk_space) + " " + (if use_ssd then "SSD" else "HDD")
    docker: "${docker}"
    cpu: "${cpu}"
    memory: "${memory}"
  }

}

