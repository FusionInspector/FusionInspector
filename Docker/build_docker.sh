#!/bin/bash

set -ev

VERSION=`cat VERSION.txt`

docker build --progress=plain -t trinityctat/fusioninspector:${VERSION} .
docker build -t trinityctat/fusioninspector:latest .

