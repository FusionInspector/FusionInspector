#!/bin/bash

set -ev

VERSION=`cat VERSION.txt`

docker push trinityctat/fusioninspector:${VERSION}
docker push trinityctat/fusioninspector:latest

