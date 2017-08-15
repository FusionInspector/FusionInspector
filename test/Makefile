

left_fq="test.reads_1.fastq.gz"
right_fq="test.reads_2.fastq.gz"


test:
	./runMe.pl --out_dir FusionInspector-outdir --left_fq ${left_fq} --right_fq ${right_fq}


all: test fusion_contigs_only incl_Trinity annotate coding_effect single_reads single_reads_fc_only


fusion_contigs_only:
	./runMe.pl --out_dir FusionInspector-fusionContigOnly --fusion_contigs_only --left_fq ${left_fq} --right_fq ${right_fq}

incl_Trinity:
	./runMe.pl --out_dir FusionInspector-fusionContigOnly-incl_Trinity --fusion_contigs_only --include_Trinity --left_fq ${left_fq} --right_fq ${right_fq}

annotate: fusion_contigs_only
	./runMe.pl --out_dir FusionInspector-fusionContigOnly --fusion_contigs_only --annotate --left_fq ${left_fq} --right_fq ${right_fq}

coding_effect: fusion_contigs_only
	./runMe.pl --out_dir FusionInspector-fusionContigOnly --fusion_contigs_only --examine_coding_effect --left_fq ${left_fq} --right_fq ${right_fq}


single_reads:
	./runMe.pl --out_dir FusionInspector-single_reads --left_fq ${left_fq}


single_reads_fc_only:
	./runMe.pl --out_dir FusionInspector-single_reads-fc_only --left_fq ${left_fq} --fusion_contigs_only


clean:
	rm -rf ./FusionInspector-fusionContigOnly ./FusionInspector-fusionContigOnly-incl_Trinity
	rm -rf ./FusionInspector-outdir
	rm -rf ./FusionInspector-single_reads ./FusionInspector-single_reads-fc_only


