

left_fq="test.reads_1.fastq.gz"
right_fq="test.reads_2.fastq.gz"


test:
	./runMe.pl --output_dir FusionInspector-outdir --left_fq ${left_fq} --right_fq ${right_fq} --vis


all: test test_incl_brkpt_plots test_incl_cosmic_like_preds fusion_contigs_only incl_Trinity annotate coding_effect samples_file_PE spanning_only_no_junction no_shrink_introns multiple_read_sets single_reads single_reads_fc_only



test_incl_brkpt_plots:
	./runMe.pl --output_dir FusionInspector-outdir --left_fq ${left_fq} --right_fq ${right_fq} --vis --incl_microH_expr_brkpt_plots

test_incl_cosmic_like_preds:
	./runMe.pl --output_dir FusionInspector-outdir --left_fq ${left_fq} --right_fq ${right_fq} --vis --predict_cosmic_like


test_incl_STAR_xtra_params:
	./runMe.pl --output_dir FusionInspector-STARxtraparams --left_fq ${left_fq} --right_fq ${right_fq} --STAR_xtra_params \"--limitSjdbInsertNsj 1005121\"


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

spanning_only_no_junction:
	rm -rf FusionInspector-spanningOnlyNoJunction data/spanning_only_no_junction
	mkdir -p data/spanning_only_no_junction
	printf "%s\n" \
		"&BT474f@HWI-EAS418:6:68:486:1723" \
		"&BT474f@HWI-EAS418:3:100:772:1835" \
		"&BT474f@HWI-EAS418:3:72:1391:1561" \
		"&BT474f@HWI-EAS418:5:97:422:825" \
		"&BT474f@HWI-EAS418:5:12:1606:926" \
		"&BT474f@HWI-EAS418:3:55:232:441" \
		> data/spanning_only_no_junction/spanning_only.read_names
	perl extract_fusion_supporting_reads.pl --left data/BT474f_1.fq.gz --right data/BT474f_2.fq.gz --fusion CTD-2328D6.1--MT-CO1 --support_type S --read_names_file data/spanning_only_no_junction/spanning_only.read_names --out_prefix data/spanning_only_no_junction/CTD-2328D6.1--MT-CO1.spanning_only
	printf "CTD-2328D6.1--MT-CO1\n" > data/spanning_only_no_junction/fusions.txt
	../FusionInspector --fusions data/spanning_only_no_junction/fusions.txt --genome_lib ${CTAT_GENOME_LIB} --out_prefix finspector --output_dir FusionInspector-spanningOnlyNoJunction --left_fq data/spanning_only_no_junction/CTD-2328D6.1--MT-CO1.spanning_only_1.fq --right_fq data/spanning_only_no_junction/CTD-2328D6.1--MT-CO1.spanning_only_2.fq --fusion_contigs_only --vis --max_sensitivity --CPU 4
	python3 validate_spanning_only_no_junction.py FusionInspector-spanningOnlyNoJunction CTD-2328D6.1--MT-CO1


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
	rm -rf ./FusionInspector-spanningOnlyNoJunction ./data/spanning_only_no_junction
	rm -rf ./FusionInspector-no_shrink_introns
	rm -rf ./FusionInspector-by-docker*
	rm -rf ./FusionInspector-multreadsets-outdir
