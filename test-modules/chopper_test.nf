#!/usr/bin/env nextflow

include { CHOPPER } from "../modules/nf-core/chopper/main.nf"
include { convert_to_nfcore_tuple } from "../subworkflows/local/utils/nf-core-compatibility.nf"
// note: path of module is relative to the directory containing this file! (./testing/)

// -------------------------
// Parameters
// -------------------------
// note: paths are relative to the workflow directory (from where nextflow is run)
params.reads = "${projectDir}/../assets/examples/data/mock9_sub10k.fastq.gz"

// -------------------------
// Workflow
// -------------------------
workflow {
    // testing convert_to_nfcore_tuple
    ch_reads = convert_to_nfcore_tuple(params.reads)

    contam_ref = Channel.of([])

    CHOPPER(ch_reads, contam_ref)
}
    