#!/usr/bin/env nextflow

include { convert_to_nfcore_tuple } from '../subworkflows/local/utils/nf-core-compatibility.nf'
include { EMU_DOWNLOADDB } from "../modules/local/emu/downloaddb/main.nf"
include { EMU_ABUNDANCE } from "../modules/local/emu/main.nf"

// -------------------------
// Parameters
// -------------------------
params.reads = "${projectDir}/../assets/examples/other_tools_files/emu_full_length.fa"
// params.emu_db = "${projectDir}/../databases/emu"

// -------------------------
// Workflow
// -------------------------
workflow {
    
    // reads_ch = convert_to_nfcore_tuple(params.reads)

    // output_dir = Channel.fromPath(params.output_dir)
    
    EMU_DOWNLOADDB() // test download db

    // EMU_ABUNDANCE(reads_ch)  
    // ch_versions = ch_versions.mix(EMU_ABUNDANCE.out.versions.first())
}
