#!/usr/bin/env nextflow

// run with: nextflow run test-modules/metagenomescope_test.nf
include { convert_to_nfcore_tuple } from '../subworkflows/local/utils/nf-core-compatibility.nf'
include { METAGENOMESCOPE } from "../modules/local/metagenomescope/main.nf"
// note: path of module is relative to the directory containing this file! (./testing/)

// -------------------------
// Parameters
// -------------------------
// note: paths are relative to the directory containing this file
params.graph = "${projectDir}/../assets/data/other_tools_files/more_files/agb/flye_yeast.gv"
params.info = "${projectDir}/../assets/data/other_tools_files/more_files/agb/flye_yeast_assembly_info.txt"

// -------------------------
// Workflow
// -------------------------
workflow {

    graph = convert_to_nfcore_tuple(params.graph)
    info = convert_to_nfcore_tuple(params.info) // pass a file(null) if no info file is provided
    
    METAGENOMESCOPE(graph, info)
}
    