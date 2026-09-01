#!/usr/bin/env bash
set -euo pipefail

echo "Pulling container images for somatic variant calling pipeline..."
echo ""

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
IMAGES_DIR="$SCRIPT_DIR/../resources/images"
RUNTIME="${1:-}"

mkdir -p "$IMAGES_DIR"

if [ -z "$RUNTIME" ]; then
    if command -v apptainer >/dev/null 2>&1; then
        RUNTIME="apptainer"
    elif command -v singularity >/dev/null 2>&1; then
        RUNTIME="singularity"
    elif command -v docker >/dev/null 2>&1; then
        RUNTIME="docker"
    else
        echo "ERROR: No container runtime found (apptainer, singularity, or docker required)" >&2
        exit 1
    fi
fi

case "$RUNTIME" in
    apptainer|singularity)
        if ! command -v "$RUNTIME" >/dev/null 2>&1; then
            echo "ERROR: $RUNTIME is not available on PATH" >&2
            exit 1
        fi

        echo "Runtime: $RUNTIME"
        echo "Images will be stored in: $IMAGES_DIR"
        echo ""

        declare -A SIF_IMAGES=(
            ["docker://broadinstitute/gatk:4.4.0.0"]="gatk.sif"
            ["docker://quay.io/biocontainers/strelka:2.9.10--0"]="strelka.sif"
            ["docker://ensemblorg/ensembl-vep@sha256:f354dd8d09073e4d943acbbd02f5eb234a9d9e9d444371c1c349910f2123de11"]="ensembl-vep-116.0.sif"
            ["docker://brentp/somalier@sha256:a99d59a80bb24d2d9d4bffc36d891cb6c935fb3ecfb61c240e97a20e19ea7916"]="somalier-0.3.3.sif"
        )

        for docker_uri in "${!SIF_IMAGES[@]}"; do
            local_name="${SIF_IMAGES[$docker_uri]}"
            output="$IMAGES_DIR/$local_name"
            if [ -f "$output" ]; then
                echo "Image already exists: $local_name"
            else
                echo "Pulling $docker_uri -> $local_name"
                "$RUNTIME" pull "$output" "$docker_uri"
            fi
        done
        ;;
    docker)
        if ! command -v docker >/dev/null 2>&1; then
            echo "ERROR: docker is not available on PATH" >&2
            exit 1
        fi

        echo "Runtime: docker"
        echo ""

        DOCKER_IMAGES=(
            "broadinstitute/gatk:4.4.0.0"
            "quay.io/biocontainers/strelka:2.9.10--0"
            "ensemblorg/ensembl-vep@sha256:f354dd8d09073e4d943acbbd02f5eb234a9d9e9d444371c1c349910f2123de11"
            "brentp/somalier@sha256:a99d59a80bb24d2d9d4bffc36d891cb6c935fb3ecfb61c240e97a20e19ea7916"
        )

        for image in "${DOCKER_IMAGES[@]}"; do
            echo "Pulling $image"
            docker pull "$image"
        done
        ;;
    *)
        echo "ERROR: Unknown runtime '$RUNTIME' (expected apptainer, singularity, or docker)" >&2
        exit 1
        ;;
esac

echo ""
echo "Done!"
