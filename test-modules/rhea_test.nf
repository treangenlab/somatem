#!/usr/bin/env nextflow

include { RHEA } from "../modules/local/rhea/main.nf"
// note: path of module is relative to the directory containing this file! (./testing/)

// -------------------------
// Parameters
// -------------------------
// note: paths are relative to the directory containing this file
params.reads = "${projectDir}/../assets/examples/other_tools_files/rhea/*.fasta"

// -------------------------
// Workflow
// -------------------------
workflow {

    reads_ch = Channel.fromPath(params.reads).collect()
    RHEA(reads_ch)
}
    