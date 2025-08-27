#!/usr/bin/env nextflow

// include { convert_to_nfcore_tuple } from '../subworkflows/local/utils/nf-core-compatibility.nf'
include { AGB } from "../modules/local/agb/main.nf"
// note: path of module is relative to the directory containing this file! (./testing/)

// -------------------------
// Parameters
// -------------------------
// note: paths are relative to the directory containing this file
params.assembly_dir = "${projectDir}/../archive/rhea_results/metaflye/"

// -------------------------
// Workflow
// -------------------------
workflow {

    assembly_dir = Channel.fromPath(params.assembly_dir)
    AGB(assembly_dir)
}
    