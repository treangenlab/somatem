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

    ch_versions = Channel.empty() // collect versions from all modules
    
    RawNanoPlot(reads_ch) // initial QC
    ch_versions = ch_versions.mix(RawNanoPlot.out.versions.first())
    
    runHostile(reads_ch, params.hostile_index) // host contamination removal // TODO: make conditional param to enable/disable
    ch_versions = ch_versions.mix(runHostile.out.versions)

    CHOPPER(runHostile.out.dehosted_reads, contam_ref) // quality filtering; future contam removal)
    ch_versions = ch_versions.mix(CHOPPER.out.versions.first())
    
    FinalNanoPlot(CHOPPER.out.fastq) // final QC

    emit:
    clean_reads = CHOPPER.out.fastq
    versions = ch_versions
}
