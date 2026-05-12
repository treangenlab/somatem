#!/usr/bin/env nextflow

include { convert_to_nfcore_tuple } from '../subworkflows/local/utils/nf-core-compatibility.nf'
include { RHEA } from "../modules/local/rhea/main.nf"
// note: path of module is relative to the directory containing this file! (./testing/)

// -------------------------
// Parameters
// -------------------------
// note: paths are relative to the directory containing this file
params.reads = "${projectDir}/../assets/data/other_tools_files/rhea/*.fasta"

// -------------------------
// Workflow
// -------------------------
workflow {

    // reads_ch = channel.fromPath(params.reads).collect()
    reads_ch = convert_to_nfcore_tuple(params.reads)

    RHEA(reads_ch)
}
    