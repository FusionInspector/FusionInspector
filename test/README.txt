To run the example you can execute the following command:


  ../FusionInspector --fusions test_fusions.list,test_fusions.list2,test_fusions.list3 \
                     --genome_lib $(/PATH/TO/CTAT_GENOME_LIB) \
                     --left_fq reads.left.simPE.fq.gz \
                     --right reads.right.simPE.fq.gz \
                     --out_dir Fusion_Inspector-STAR \
                     --out_prefix finspector \
                     --vis


Note, you must have the CTAT_GENOME_LIB installed as per: https://github.com/FusionInspector/FusionInspector/wiki/installing-FusionInspector#data-requirements


