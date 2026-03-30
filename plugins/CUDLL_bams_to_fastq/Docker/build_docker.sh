#!/bin/bash

set -ev

VERSION=`cat VERSION.txt`

docker build -t trinityctat/cudll-to-fastq:latest .
docker build -t trinityctat/cudll-to-fastq:${VERSION} .

