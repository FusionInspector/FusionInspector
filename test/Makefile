

left_fq="test.reads_1.fastq.gz"
right_fq="test.reads_2.fastq.gz"


test:
	./runMe.pl --out_dir FusionInspector-outdir --left_fq ${left_fq} --right_fq ${right_fq}


all: test fusion_contigs_only incl_Trinity annotate coding_effect single_reads single_reads_fc_only samples_file_PE


fusion_contigs_only:
	./runMe.pl --out_dir FusionInspector-fusionContigOnly --fusion_contigs_only --left_fq ${left_fq} --right_fq ${right_fq} --extract_fusion_reads_file FusionInspector-fusionContigOnly/fusion_reads --vis

incl_Trinity:
	./runMe.pl --out_dir FusionInspector-fusionContigOnly-incl_Trinity --fusion_contigs_only --include_Trinity --left_fq ${left_fq} --right_fq ${right_fq} --vis --write_intermediate_results

annotate: fusion_contigs_only
	./runMe.pl --out_dir FusionInspector-fusionContigOnly --fusion_contigs_only --annotate --left_fq ${left_fq} --right_fq ${right_fq}

coding_effect: fusion_contigs_only
	./runMe.pl --out_dir FusionInspector-fusionContigOnly --fusion_contigs_only --examine_coding_effect --left_fq ${left_fq} --right_fq ${right_fq}


single_reads:
	./runMe.pl --out_dir FusionInspector-single_reads --left_fq ${left_fq} --extract_fusion_reads_file FusionInspector-single_reads/fusion_reads


single_reads_fc_only:
	./runMe.pl --out_dir FusionInspector-single_reads-fc_only --left_fq ${left_fq} --fusion_contigs_only --extract_fusion_reads_file FusionInspector-single_reads/fusion_reads 


samples_file_PE:
	./runMe.pl --out_dir FusionInspector-pe_samples --samples_file data/samples.txt --fusion_contigs_only --extract_fusion_reads_file FusionInspector-pe_samples/fusion_reads


clean:
	rm -rf ./FusionInspector-fusionContigOnly ./FusionInspector-fusionContigOnly-incl_Trinity
	rm -rf ./FusionInspector-outdir
	rm -rf ./FusionInspector-single_reads ./FusionInspector-single_reads-fc_only
	rm -rf ./FusionInspector-pe_samples
	rm -rf ./FusionInspector-by-docker

