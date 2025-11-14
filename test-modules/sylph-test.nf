#!/usr/bin/env nextflow

include { SYLPH } from "../modules/local/sylph_profile/main.nf"

// -------------------------
// Parameters
// -------------------------
params.reads = './assets/examples/sylph/o157_reads.fastq.gz'
params.reference = './assets/examples/sylph/*.fasta.gz' // tiny : 3 genomes only
// params.reference = './databases/v0.3-c1000-gtdb-r214.syldb' // small database
params.threads = 4

// -------------------------
// Workflow
// -------------------------
workflow {

    reads = Channel.fromPath(params.reads)
    reference = Channel.fromPath(params.reference).collect()

    SYLPH(reads, reference, params.threads)
    // reads.view {read -> "read: ${read.baseName}"}
    reference.view {ref -> "reference: ${ref.baseName}"}
}
    