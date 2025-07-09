#!/usr/bin/env nextflow

include { CENTRIFUGER_CLASSIFY } from "../modules/local/centrifuger/main.nf"
// note: path of module is relative to the directory containing this file! (./testing/)

// -------------------------
// Parameters
// -------------------------
// note: paths are relative to the workflow directory (from where nextflow is run)
params.reads = "${projectDir}/../examples/data/46_1_sub10k.fastq.gz"
params.index_prefix = "${projectDir}/../databases/legionella_idx/cfr_ref_idx"
params.threads = 4

// -------------------------
// Workflow
// -------------------------
workflow {

    reads = Channel.fromPath(params.reads)
    index_prefix = Channel.fromPath(params.index_prefix)
    threads = Channel.of(params.threads)

    CENTRIFUGER_CLASSIFY(reads, index_prefix, threads)
}
    