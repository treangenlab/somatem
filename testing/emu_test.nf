#!/usr/bin/env nextflow


include { EMU_ABUNDANCE } from "../modules/local/emu/main.nf"

// -------------------------
// Parameters
// -------------------------
params.reads = './examples/data/emu_full_length.fa'

// -------------------------
// Workflow
// -------------------------
workflow {
    
    reads = Channel.fromPath(params.reads)
    meta = Channel.of([id: 'test', single_end: false]) // meta is required by EMU, initialize with dummy values
    temp_meta_reads = tuple(meta, reads)

    // output_dir = Channel.fromPath(params.output_dir)
    
    // temp_meta_reads.view { meta, reads -> "meta: $meta, reads: $reads" }

    EMU_ABUNDANCE(temp_meta_reads)  
    // ch_versions = ch_versions.mix(EMU_ABUNDANCE.out.versions.first())
}
