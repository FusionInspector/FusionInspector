#!/bin/bash

set -ve

## assumes we've run the docker test already

CTAT_GENOME_LIB="GRCh37_gencode_v19_CTAT_lib_Aug152019.plug-n-play"


VERSION=`cat VERSION.txt`

# run FusionInspector

TESTDIR=test
fusion_files_list="${TESTDIR}/fusion_targets.A.txt,${TESTDIR}/fusion_targets.B.txt,${TESTDIR}/fusion_targets.C.txt"
left_fq="${TESTDIR}/test.reads_1.fastq.gz"
right_fq="${TESTDIR}/test.reads_2.fastq.gz"

cd ../ && singularity exec -e Docker/FusionInspector.v${VERSION}.simg FusionInspector  \
       --fusions $fusion_files_list \
       --output_dir ${TESTDIR}/FusionInspector-by-singularity \
       --left_fq ${left_fq} \
       --right_fq ${right_fq} \
       --genome_lib ${CTAT_GENOME_LIB}/ctat_genome_lib_build_dir  \
       --out_prefix finspector \
       --vis \
       --include_Trinity \
       --examine_coding_effect \
       --extract_fusion_reads_file ${TESTDIR}/FusionInspector-by-singularity/fusion_reads 

