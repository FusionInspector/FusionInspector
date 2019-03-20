#!/bin/bash

set -ev

VERSION=`cat VERSION.txt`

docker build -t trinityctat/fusioninspector:${VERSION} .
docker build -t trinityctat/fusioninspector:latest .

