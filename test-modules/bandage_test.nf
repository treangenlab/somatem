#!/usr/bin/env nextflow

include { convert_to_nfcore_tuple } from '../subworkflows/local/utils/nf-core-compatibility.nf'
include { BANDAGENG_IMAGE } from "../modules/local/bandageng/image/main.nf"
// note: path of module is relative to the directory containing this file! (./testing/)

// -------------------------
// Parameters
// -------------------------
// note: paths are relative to the directory containing this file
params.graph_file = "${projectDir}/../assets/data/intermediate_files/rhea-t0-t1/metaflye/assembly_graph.gfa"
params.colour_csv = "${projectDir}/../assets/data/intermediate_files/rhea-t0-t1/Bandage_metadata.csv"

// -------------------------
// Workflow
// -------------------------
workflow {

    graph_file = convert_to_nfcore_tuple(params.graph_file)
    colour_csv = convert_to_nfcore_tuple(params.colour_csv)

    BANDAGENG_IMAGE(graph_file, colour_csv)
}
    