#!/usr/bin/env nextflow

include { TOOL } from "../modules/local/module_template.nf"
// note: path of module is relative to the directory containing this file! (./testing/)

// -------------------------
// Parameters
// -------------------------
// note: paths are relative to the directory containing this file
params.reads = "${projectDir}/../examples/data/emu_full_length.fa"
params.database = "${projectDir}/../examples/data/lemur_report.csv"
params.taxonomy = "${projectDir}/../examples/data/lemur_report.csv"
params.rank = 'species'

// -------------------------
// Workflow
// -------------------------
workflow {

    reads = Channel.fromPath(params.reads)
    database = Channel.fromPath(params.database)
    taxonomy = Channel.fromPath(params.taxonomy)
    rank = Channel.of(params.rank)

    TOOL(reads, database, taxonomy, rank)
}
    