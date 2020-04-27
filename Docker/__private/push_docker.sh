#!/bin/bash

set -ev

VERSION=`cat VERSION.txt`

docker push bhaastestdockers/fusioninspector:${VERSION}


