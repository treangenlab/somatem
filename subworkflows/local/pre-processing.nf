#!/usr/bin/env nextflow

include { convert_to_nfcore_tuple } from './utils/nf-core-compatibility.nf'

include { NANOPLOT as RawNanoPlot; NANOPLOT as FinalNanoPlot } from '../../modules/nf-core/nanoplot/main.nf'
include { runHostile } from './runHostile.nf'
include { CHOPPER } from '../../modules/nf-core/chopper/main.nf'


// -------------------------
// Workflow Definition
// -------------------------

workflow PREPROCESSING {

    take:
    reads_ch
    contam_ref

    main:
    RawNanoPlot(reads_ch) // initial QC
    
    runHostile(reads_ch, params.hostile_index) // host contamination removal // TODO: make conditional param to enable/disable
    CHOPPER(runHostile.out, contam_ref) // quality filtering; future contam removal)
    
    FinalNanoPlot(CHOPPER.out.fastq) // final QC

    emit:
    clean_reads = CHOPPER.out.fastq
}
