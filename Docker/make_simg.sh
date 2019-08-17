#!/bin/bash

VERSION=`cat VERSION.txt`

singularity build FusionInspector.v${VERSION}.simg docker://trinityctat/fusioninspector:$VERSION

singularity exec -e FusionInspector.v${VERSION}.simg FusionInspector --version


