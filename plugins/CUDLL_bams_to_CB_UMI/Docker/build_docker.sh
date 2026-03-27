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

docker build \
    -t "${IMAGE_NAME}:latest" \
    -t "${IMAGE_NAME}:${IMAGE_VERSION}" \
    "${SCRIPT_DIR}"

echo ""
echo "Docker image built successfully: ${IMAGE_NAME}:latest"
echo "Docker image built successfully: ${IMAGE_NAME}:${IMAGE_VERSION}"
