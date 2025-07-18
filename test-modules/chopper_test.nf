#!/usr/bin/env nextflow

include { CHOPPER } from "../modules/nf-core/chopper/main.nf"
// note: path of module is relative to the directory containing this file! (./testing/)

// -------------------------
// Parameters
// -------------------------
// note: paths are relative to the workflow directory (from where nextflow is run)
params.reads = "${projectDir}/../examples/data/46_1_sub10k.fastq.gz"

// -------------------------
// Workflow
// -------------------------
workflow {
    // testing
    ch_reads = Channel.fromPath(params.reads)
                    .map { r ->
                        def meta = [:] // Use dummy values; meta is required by nf-core modules
                        meta.id = "test"
                        meta.single_end = false
                        return [meta, r] }

    contam_ref = Channel.of([])

    ch_reads.view { r -> "tuple: ${r.class}, meta: ${r[0].class}, reads: ${r[1].class}" }
    contam_ref.view { r -> "contam_ref: $r.class" }

    CHOPPER(ch_reads, contam_ref)
}
    