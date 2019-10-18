#!/bin/bash

set -ve

if [ -z ${CTAT_GENOME_LIB} ]; then
    echo "Error, must have CTAT_GENOME_LIB env var set"
    exit 1
fi


VERSION=`cat VERSION.txt`

# run FusionInspector

TESTDIR=test
fusion_files_list="${TESTDIR}/fusion_targets.A.txt,${TESTDIR}/fusion_targets.B.txt,${TESTDIR}/fusion_targets.C.txt"
left_fq="${TESTDIR}/test.reads_1.fastq.gz"
right_fq="${TESTDIR}/test.reads_2.fastq.gz"

cd ../ && singularity exec -e -B ${CTAT_GENOME_LIB}:/ctat_genome_lib:ro \
          Docker/FusionInspector.v${VERSION}.simg FusionInspector  \
       --fusions $fusion_files_list \
       --output_dir ${TESTDIR}/FusionInspector-by-singularity \
       --left_fq ${left_fq} \
       --right_fq ${right_fq} \
       --genome_lib /ctat_genome_lib  \
       --out_prefix finspector \
       --vis \
       --include_Trinity \
       --examine_coding_effect \
       --extract_fusion_reads_file ${TESTDIR}/FusionInspector-by-singularity/fusion_reads 

