#!/bin/bash

set -ve

if [ -z $CTAT_GENOME_LIB ]; then
    echo "Error, must have CTAT_GENOME_LIB env var set"
    exit 1
fi


VERSION=`cat VERSION.txt`

# run FusionInspector

TESTDIR=/data/test
fusion_files_list="${TESTDIR}/fusion_targets.A.txt,${TESTDIR}/fusion_targets.B.txt,${TESTDIR}/fusion_targets.C.txt"
left_fq="${TESTDIR}/test.reads_1.fastq.gz"
right_fq="${TESTDIR}/test.reads_2.fastq.gz"

docker run -v `pwd`/../:/data -v ${CTAT_GENOME_LIB}:/ctat_genome_lib \
       --rm trinityctat/fusioninspector:${VERSION} FusionInspector  \
       --fusions $fusion_files_list \
       --output_dir ${TESTDIR}/FusionInspector-by-docker \
       --left_fq ${left_fq} \
       --right_fq ${right_fq} \
       --genome_lib /ctat_genome_lib  \
       --out_prefix finspector \
       --vis \
       --include_Trinity \
       --examine_coding_effect \
       --extract_fusion_reads_file ${TESTDIR}/FusionInspector-by-docker/fusion_reads 

