#!/usr/bin/env nextflow

include { GDOWN } from '../../modules/local/gdown/main.nf'

workflow {
    GDOWN()
}
    