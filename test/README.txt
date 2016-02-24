To run the example you can execute the following command:


  ../FusionInspector --fusions test_fusions.list,test_fusions.list2,test_fusions.list3 \
                     --genome_lib $(/PATH/TO/CTAT_GENOME_LIB) \
                     --left_fq reads.left.simPE.fq.gz \
                     --right reads.right.simPE.fq.gz \
                     --out_dir Fusion_Inspector-STAR \
                     --out_prefix finspector \
                     --align_utils STAR --prep_for_IGV --no_cleanup


Note, you must have the CTAT_GENOME_LIB installed as per: https://github.com/FusionFilter/FusionFilter/wiki


Alternatively to running the above command, you can 

   export CTAT_GENOME_LIB=$(/PATH/TO/CTAT_GENOME_LIB)

   ./runMe.pl


and it will execute automatically for you, but this does require that you have the env var set as above.

