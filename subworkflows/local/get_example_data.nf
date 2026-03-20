#!/usr/bin/env nextflow

include { GDOWN } from '../../modules/local/gdown/main.nf'

params.outdir = params.outdir ?: "${launchDir}/somatem_example_data"

workflow {
    GDOWN()
}