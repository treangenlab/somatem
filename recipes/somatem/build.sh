#!/usr/bin/env bash
set -euo pipefail

# Get the path of the conda package within the conda environment
# $PREFIX = path of the conda environment: conda install --prefix <path> somatem
PIPELINE_DIR="${PREFIX}/share/${PKG_NAME}-${PKG_VERSION}"

mkdir -p "${PREFIX}/bin"
mkdir -p "${PIPELINE_DIR}"

# copy pipeline executable
cp bin/somatem "${PREFIX}/bin/${PKG_NAME}"
chmod +x "${PREFIX}/bin/${PKG_NAME}"

# substitute placeholder with actual pipeline directory in the executable file
sed -i "s|@@PIPELINE_DIR@@|${PIPELINE_DIR}|g" "${PREFIX}/bin/${PKG_NAME}"

# copy pipeline directories
cp -r assets "${PIPELINE_DIR}/"
cp -r bin "${PIPELINE_DIR}/"
cp -r conf "${PIPELINE_DIR}/"
cp -r docs "${PIPELINE_DIR}/"
cp -r modules "${PIPELINE_DIR}/"
cp -r subworkflows "${PIPELINE_DIR}/"
cp -r workflows "${PIPELINE_DIR}/"
cp -r tests "${PIPELINE_DIR}/"

# copy pipeline files
cp CITATIONS.md "${PIPELINE_DIR}/"
cp LICENSE "${PIPELINE_DIR}/"
cp README.md "${PIPELINE_DIR}/"
cp main.nf "${PIPELINE_DIR}/"
cp modules.json "${PIPELINE_DIR}/"
cp nextflow.config "${PIPELINE_DIR}/"
cp nf-test.config "${PIPELINE_DIR}/"