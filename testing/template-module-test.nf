#!/usr/bin/env nextflow

include { TOOL } from "../modules/local/module_template.nf"

// -------------------------
// Parameters
// -------------------------
params.reads = './examples/data/emu_full_length.fa'
params.database = './examples/data/lemur_report.csv'
params.other = './examples/data/lemur_report.csv'

// -------------------------
// Workflow
// -------------------------
workflow {
    TOOL(params.reads, params.database, params.other)
}
    