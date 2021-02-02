

left_fq="test.reads_1.fastq.gz"
right_fq="test.reads_2.fastq.gz"


test:
	./runMe.pl --output_dir FusionInspector-outdir --left_fq ${left_fq} --right_fq ${right_fq} --vis


all: test fusion_contigs_only incl_Trinity annotate coding_effect samples_file_PE no_shrink_introns multiple_read_sets single_reads single_reads_fc_only


no_shrink_introns:
	./runMe.pl --output_dir FusionInspector-no_shrink_introns --left_fq ${left_fq} --right_fq ${right_fq} --no_shrink_introns 

fusion_contigs_only:
	./runMe.pl --output_dir FusionInspector-fusionContigOnly --fusion_contigs_only --left_fq ${left_fq} --right_fq ${right_fq} --extract_fusion_reads_file FusionInspector-fusionContigOnly/fusion_reads --vis

fusion_contigs_only_skip_EM:
	./runMe.pl --output_dir FusionInspector-fusionContigOnly-skip_EM --fusion_contigs_only --left_fq ${left_fq} --right_fq ${right_fq} --extract_fusion_reads_file FusionInspector-fusionContigOnly-skip_EM/fusion_reads --skip_EM


multiple_read_sets:
	./runMe.pl --output_dir FusionInspector-multreadsets-outdir --left_fq test.readsA_1.fastq.gz,test.readsB_1.fastq.gz --right_fq test.readsA_2.fastq.gz,test.readsB_2.fastq.gz


incl_Trinity:
	./runMe.pl --output_dir FusionInspector-fusionContigOnly-incl_Trinity --fusion_contigs_only --include_Trinity --left_fq ${left_fq} --right_fq ${right_fq} --vis --write_intermediate_results

annotate: fusion_contigs_only
	./runMe.pl --output_dir FusionInspector-fusionContigOnly --fusion_contigs_only --annotate --left_fq ${left_fq} --right_fq ${right_fq}

coding_effect: fusion_contigs_only
	./runMe.pl --output_dir FusionInspector-fusionContigOnly --fusion_contigs_only --examine_coding_effect --left_fq ${left_fq} --right_fq ${right_fq}


single_reads:
	./runMe.pl --output_dir FusionInspector-single_reads --left_fq ${left_fq} --extract_fusion_reads_file FusionInspector-single_reads/fusion_reads


single_reads_fc_only:
	./runMe.pl --output_dir FusionInspector-single_reads-fc_only --left_fq ${left_fq} --fusion_contigs_only --extract_fusion_reads_file FusionInspector-single_reads/fusion_reads 


samples_file_PE:
	./runMe.pl --output_dir FusionInspector-pe_samples --samples_file data/samples.txt --fusion_contigs_only --extract_fusion_reads_file FusionInspector-pe_samples/fusion_reads


######################
## single cell support

single_cell_PE:
	../util/sc/prep_distributed_jobs.py --sample_sheet data/samples.txt --cells_per_job 1 --output_dir FI_SingleCellPE # cells_per_job set to 1 just for this smallt est. Default is actually 24
	../util/sc/write_sc_FI_cmds.py --batches_list_file FI_SingleCellPE.batches.list --genome_lib_dir ${CTAT_GENOME_LIB} --fusions fusion_targets.A.txt > SingleCellPE.batches.FI.cmds
	../util/sc/run_distributed_jobs_locally.py --cmds_file SingleCellPE.batches.FI.cmds --num_parallel_exec 4 --genome_lib_dir ${CTAT_GENOME_LIB}
	../util/sc/aggregate_and_deconvolve_fusion_outputs.py --batches_list_file FI_SingleCellPE.batches.list --output_prefix FI_SingleCellPE



clean:
	rm -rf ./FusionInspector-fusionContigOnly ./FusionInspector-fusionContigOnly-incl_Trinity ./FusionInspector-by-singularity
	rm -rf ./FusionInspector-outdir
	rm -rf ./FusionInspector-single_reads ./FusionInspector-single_reads-fc_only
	rm -rf ./FusionInspector-pe_samples
	rm -rf ./FusionInspector-no_shrink_introns
	rm -rf ./FusionInspector-by-docker
	rm -rf ./FusionInspector-multreadsets-outdir


