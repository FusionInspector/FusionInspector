#!/bin/bash

set -ev

VERSION=`cat VERSION.txt`

docker push trinityctat/cudll-to-fastq:latest
docker push trinityctat/cudll-to-fastq:${VERSION}

