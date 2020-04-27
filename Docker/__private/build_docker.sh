#!/bin/bash

set -ev

VERSION=`cat VERSION.txt`

docker build -t bhaastestdockers/fusioninspector:${VERSION} .


