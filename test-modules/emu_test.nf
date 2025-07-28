#!/usr/bin/env nextflow


include { EMU_ABUNDANCE } from "../modules/local/emu/main.nf"

// -------------------------
// Parameters
// -------------------------
params.reads = "${projectDir}/../examples/data/emu_full_length.fa"
params.db = "${projectDir}/../databases/emu"

// -------------------------
// Workflow
// -------------------------
workflow {
    
    reads = Channel.fromPath(params.reads)
    meta = Channel.of([id: 'test', single_end: false]) // meta is required by EMU, initialize with dummy values
    db = Channel.fromPath(params.db)
    
    // output_dir = Channel.fromPath(params.output_dir)
    
    meta.view { m -> "meta: $m" }
    reads.view { r -> "reads: $r" }

    EMU_ABUNDANCE(meta, reads, db)  
    // ch_versions = ch_versions.mix(EMU_ABUNDANCE.out.versions.first())
}
