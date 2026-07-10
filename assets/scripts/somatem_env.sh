# !/bin/bash
# This file is meant to be sourced by users to set environment variables for the Somatem pipeline 
# to change default locations for databases and conda cache
# It is not required to run the pipeline, but can be used to customize the environment.

# Set the base directory for Somatem databases. default: ~/somatem_databases
export SOMATEM_DB_DIR=/home/dbs

# Keep Nextflow environments and micromamba package downloads in the active
# Somatem installation. Set SOMATEM_HOME before sourcing to override this root.
if [[ -n "${CONDA_PREFIX:-}" ]]; then
    export SOMATEM_HOME="${SOMATEM_HOME:-${CONDA_PREFIX}/share/somatem}"
else
    export SOMATEM_HOME="${SOMATEM_HOME:-${PWD}/.somatem}"
fi
export NXF_CONDA_CACHEDIR="${NXF_CONDA_CACHEDIR:-${SOMATEM_CONDA_CACHE:-${SOMATEM_HOME}/conda_cache}}"
export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-${SOMATEM_HOME}/micromamba}"

# path to Eddy's unified databases (Ensemble analysis: species detection)
export SOMATEM_UNIFIED_DB_DIR=/home/Users/pacbio_bakeoff/data/ref_db/refseq03032025
