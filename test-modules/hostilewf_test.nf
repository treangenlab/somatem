#!/usr/bin/env nextflow

include { runHostile } from "../subworkflows/local/runHostile.nf"
include { convert_to_nfcore_tuple } from "../subworkflows/local/utils/nf-core-compatibility.nf"

// -------------------------
// Parameters
// -------------------------
params.reads = "${projectDir}/../assets/examples/data/mock9_sub10k.fastq.gz"


workflow {

    reads_ch = convert_to_nfcore_tuple(params.reads)

    runHostile(reads_ch, params.hostile_index)
}
