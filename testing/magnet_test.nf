#!/usr/bin/env nextflow

include { MAGNET } from "../modules/local/magnet/main.nf"

// -------------------------
// Parameters
// -------------------------
params.reads = './examples/data/emu_full_length.fa'
params.classification = './examples/data/lemur_report.csv'
params.output = './examples/data/lemur_report.csv'

// -------------------------
// Workflow
// -------------------------
workflow {
    MAGNET(params.reads, params.classification, params.output)
}
    