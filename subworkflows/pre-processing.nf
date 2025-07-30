#!/usr/bin/env nextflow

include { convert_to_nfcore_tuple } from './utils/nf-core-compatibility.nf'

include { NANOPLOT as RawNanoPlot; NANOPLOT as FinalNanoPlot } from '../modules/nf-core/nanoplot/main.nf'
include { runHostile } from './runHostile.nf'
include { CHOPPER } from '../modules/nf-core/chopper/main.nf'




// -------------------------
// Parameters
// -------------------------

params.input_dir   = 'examples/data'
// params.output_dir  = 'results'
// params.threads     = 64
// params.maxlength   = 30000
// params.minq        = 10
// params.minlen      = 250
params.host_index  = 'human-t2t-hla-argos985-mycob140.mmi'

// -------------------------
// Workflow Definition
// -------------------------

workflow {
    in_ch = convert_to_nfcore_tuple(params.input_dir)
    contam_ref = Channel.value(null) // empty channel for now
    
    // Debug: Check what each channel contains
    // in_ch.view { ch -> "in_ch: $ch" }
    // contam_ref.view { ref -> "contam_ref: $ref" }


    // RawNanoPlot(in_ch) // initial QC
    runHostile(in_ch, params.host_index) // host contamination removal

    // Debug: Check runHostile output
    // runHostile.out.view { ch -> "runHostile.out: ${ch[1].simpleName}" }

    CHOPPER(runHostile.out, contam_ref) // quality filtering; future contam removal)
    // FinalNanoPlot(CHOPPER.out.fastq) // final QC
}
