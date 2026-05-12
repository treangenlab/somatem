#!/usr/bin/env nextflow

include { convert_to_nfcore_tuple } from '../subworkflows/local/utils/nf-core-compatibility.nf'
include { AGB } from "../modules/local/agb/main.nf"
// note: path of module is relative to the directory containing this file! (./testing/)

// -------------------------
// Parameters
// -------------------------
// note: paths are relative to the directory containing this file
// params.assembly_dir = "${projectDir}/../archive/rhea_results/metaflye/"
params.assembly_dir = "${projectDir}/../assets/data/intermediate_files/metaflye_rhea-t0-t1/"

// -------------------------
// Workflow
// -------------------------
workflow {

    assembly_dir = channel.fromPath(params.assembly_dir)
        .map { r ->
                def meta = [:] // Use dummy values; meta is required by nf-core modules
                meta.id = r.simpleName
                meta.single_end = true
                return [meta, r] }
                
    AGB(assembly_dir)
}
    