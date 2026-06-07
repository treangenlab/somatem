# !/bin/bash
# This file is meant to be sourced by users to set environment variables for the Somatem pipeline 
# to change default locations for databases and conda cache
# It is not required to run the pipeline, but can be used to customize the environment.

# Set the base directory for Somatem databases. default: ~/somatem_databases
export SOMATEM_DB_DIR=/home/dbs

# Set the conda/micromamba cache directory for Nextflow to use. default: ~/.nextflow/cache
export NXF_CONDA_CACHEDIR=~/micromamba/nextflow-envs

# path to Eddy's unified databases (Ensemble analysis: species detection)
export SOMATEM_UNIFIED_DB_DIR=/home/Users/pacbio_bakeoff/data/ref_db/refseq03032025