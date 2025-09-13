#!/usr/bin/env nextflow

include { CENTRIFUGER_CLASSIFY } from "../modules/local/centrifuger/main.nf"
// note: path of module is relative to the directory containing this file! (./testing/)

// -------------------------
// Parameters
// -------------------------
// note: paths are relative to the workflow directory (from where nextflow is run)
params.reads = "${projectDir}/../assets/examples/data/mock9_sub10k.fastq.gz"
params.db_dir = "${projectDir}/../databases/legionella_cfr_idx/"
// full path with: cfr_ref_idx.* doesn't work.
params.threads = 4

// -------------------------
// Workflow
// -------------------------
workflow {

    reads = Channel.fromPath(params.reads)
    db_dir = Channel.fromPath(params.db_dir)
    threads = Channel.of(params.threads)

    CENTRIFUGER_CLASSIFY(reads, db_dir, threads)
}
    