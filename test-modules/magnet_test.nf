#!/usr/bin/env nextflow

// test magnet module using : nextflow run test-modules/magnet_test.nf
include { convert_to_nfcore_tuple } from '../subworkflows/local/utils/nf-core-compatibility.nf'
include { MAGNET } from "../modules/local/magnet/main.nf"

// -------------------------
// Parameters
// -------------------------
params.reads = "${projectDir}/../assets/data/other_tools_files/lemur/example-data/example.fastq"
params.classification = "${projectDir}/../assets/data/other_tools_files/lemur/example-output-ref/relative_abundance.tsv"

// -------------------------
// Workflow
// -------------------------
workflow {
    reads = convert_to_nfcore_tuple(params.reads)
    classification = channel.fromPath(params.classification)

    // reads.view()
    // classification.view()

    MAGNET(reads, classification)
}
    