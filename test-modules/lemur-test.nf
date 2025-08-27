#!/usr/bin/env nextflow

// enable dsl2 syntax
nextflow.enable.dsl = 2

include { convert_to_nfcore_tuple } from '../subworkflows/local/utils/nf-core-compatibility.nf'
include { LEMUR } from "../modules/local/lemur/main.nf"

// -------------------------
// Parameters
// -------------------------
// Note: Paths are relative to the base directory of the workflow (where nextflow is run from)

params.reads = "${projectDir}/../examples/lemur/example-data/example.fastq"

// -------------------------
// Workflow
// -------------------------
workflow {
    
    reads_ch = convert_to_nfcore_tuple(params.reads)
    
    LEMUR(reads_ch)
}