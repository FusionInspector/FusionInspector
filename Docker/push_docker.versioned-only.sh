#!/bin/bash

set -ev

VERSION=`cat VERSION.txt`

docker push trinityctat/fusioninspector:${VERSION}


