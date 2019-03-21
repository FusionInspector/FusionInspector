#!/bin/bash

set -ve

CTAT_GENOME_LIB="GRCh37_v19_CTAT_lib_Feb092018"

CTAT_GENOME_LIB_URL="https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/__genome_libs_StarFv1.3/GRCh37_v19_CTAT_lib_Feb092018.plug-n-play.tar.gz"


if [ ! -s "../${CTAT_GENOME_LIB}.tar.gz" ]; then
    wget ${CTAT_GENOME_LIB_URL} -O ../${CTAT_GENOME_LIB}.tar.gz
fi


if [ ! -d "../${CTAT_GENOME_LIB}" ]; then
    tar xvf ../${CTAT_GENOME_LIB}.tar.gz -C ../.
fi


VERSION=`cat VERSION.txt`

# run FusionInspector

TESTDIR=/data/test
fusion_files_list="${TESTDIR}/fusion_targets.A.txt,${TESTDIR}/fusion_targets.B.txt,${TESTDIR}/fusion_targets.C.txt"
left_fq="${TESTDIR}/test.reads_1.fastq.gz"
right_fq="${TESTDIR}/test.reads_2.fastq.gz"

docker run -v `pwd`/../:/data --rm trinityctat/fusioninspector:${VERSION} FusionInspector  \
       --fusions $fusion_files_list \
       --out_dir ${TESTDIR}/FusionInspector-by-docker \
       --left_fq ${left_fq} \
       --right_fq ${right_fq} \
       --genome_lib /data/${CTAT_GENOME_LIB}/ctat_genome_lib_build_dir  \
       --out_prefix finspector \
       --vis \
       --include_Trinity \
       --examine_coding_effect \
       --extract_fusion_reads_file ${TESTDIR}/FusionInspector-by-docker/fusion_reads 

