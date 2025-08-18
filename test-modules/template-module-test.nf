#!/usr/bin/env nextflow

include { convert_to_nfcore_tuple } from '../subworkflows/local/utils/nf-core-compatibility.nf'
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

    reads_ch = convert_to_nfcore_tuple(params.reads)
    database_ch = Channel.fromPath(params.database)
    taxonomy_ch = Channel.fromPath(params.taxonomy)
    rank_ch = Channel.of(params.rank)

    TOOL(reads_ch, database_ch, taxonomy_ch, rank_ch)
}
    