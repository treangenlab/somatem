#!/usr/bin/env nextflow

include { NANOPLOT } from "../modules/nf-core/nanoplot/main.nf"
include { convert_to_nfcore_tuple } from "../subworkflows/utils/nf-core-compatibility.nf"
// note: path of module is relative to the directory containing this file! (./testing/)

// -------------------------
// Parameters
// -------------------------
// note: paths are relative to the workflow directory (from where nextflow is run)
params.reads = "${projectDir}/../examples/data/46_1_sub10k.fastq.gz"

// -------------------------
// Workflow
// -------------------------
workflow {

    reads = convert_to_nfcore_tuple(params.reads)

    NANOPLOT(reads)
}
    