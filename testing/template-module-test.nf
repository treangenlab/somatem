#!/usr/bin/env nextflow

include { TOOL } from "../modules/local/module_template.nf"
// note: path of module is relative to the directory containing this file! (./testing/)

// -------------------------
// Parameters
// -------------------------
// note: paths are relative to the workflow directory (from where nextflow is run)
params.reads = './examples/data/emu_full_length.fa'
params.database = './examples/data/lemur_report.csv'
params.other = './examples/data/lemur_report.csv'

// -------------------------
// Workflow
// -------------------------
workflow {

    reads = Channel.fromPath(params.reads)
    database = Channel.fromPath(params.database)
    other = Channel.fromPath(params.other)

    TOOL(reads, database, other)
}
    