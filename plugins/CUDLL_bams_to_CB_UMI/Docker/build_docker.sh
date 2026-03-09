#!/bin/bash

set -e

docker build -t trinityctat/cudll-to-cb-umi:latest .

echo ""
echo "Docker image built successfully: trinityctat/cudll-to-cb-umi:latest"
