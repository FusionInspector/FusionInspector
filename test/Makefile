
test:
	./runMe.pl --out_dir FusionInspector-outdir


all: test fusion_contigs_only incl_Trinity annotate coding_effect


fusion_contigs_only:
	./runMe.pl --out_dir FusionInspector-fusionContigOnly --fusion_contigs_only

incl_Trinity: fusion_contigs_only
	./runMe.pl --out_dir FusionInspector-fusionContigOnly --fusion_contigs_only --include_Trinity

annotate: fusion_contigs_only
	./runMe.pl --out_dir FusionInspector-fusionContigOnly --fusion_contigs_only --annotate

coding_effect: fusion_contigs_only
	./runMe.pl --out_dir FusionInspector-fusionContigOnly --fusion_contigs_only --examine_coding_effect


clean:
	rm -rf ./FusionInspector-fusionContigOnly
	rm -rf ./FusionInspector-outdir

