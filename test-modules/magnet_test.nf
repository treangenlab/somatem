#!/usr/bin/env nextflow

include { MAGNET } from "../modules/local/magnet/main.nf"

// -------------------------
// Parameters
// -------------------------
params.reads = "${projectDir}/../examples/lemur/example-data/example.fastq"
params.classification = "${projectDir}/../examples/lemur/example-output-ref/relative_abundance.tsv"
// /home/pbk1/somatem/examples/lemur/example-output-ref/relative_abundance.tsv

// -------------------------
// Workflow
// -------------------------
workflow {

    reads = Channel.fromPath(params.reads)
    classification = Channel.fromPath(params.classification)

    // reads.view()
    // classification.view()

    MAGNET(reads, classification)
}
    