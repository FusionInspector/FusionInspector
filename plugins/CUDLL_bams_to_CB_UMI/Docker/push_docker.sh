#!/bin/bash

set -euo pipefail

IMAGE_NAME="${IMAGE_NAME:-trinityctat/cudll-to-cb-umi}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VERSION_FILE="${SCRIPT_DIR}/VERSION.txt"

if [[ ! -f "${VERSION_FILE}" ]]; then
    echo "ERROR: VERSION.txt not found at ${VERSION_FILE}" >&2
    exit 1
fi

IMAGE_VERSION="$(tr -d '[:space:]' < "${VERSION_FILE}")"

if [[ -z "${IMAGE_VERSION}" ]]; then
    echo "ERROR: VERSION.txt is empty" >&2
    exit 1
fi

docker push "${IMAGE_NAME}:latest"
docker push "${IMAGE_NAME}:${IMAGE_VERSION}"

echo ""
echo "Docker image pushed successfully: ${IMAGE_NAME}:latest"
echo "Docker image pushed successfully: ${IMAGE_NAME}:${IMAGE_VERSION}"
