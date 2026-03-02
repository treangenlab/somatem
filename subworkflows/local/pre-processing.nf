#!/usr/bin/env nextflow

include { convert_to_nfcore_tuple } from './utils/nf-core-compatibility.nf'

include { NANOPLOT as RawNanoPlot; NANOPLOT as FinalNanoPlot } from '../../modules/nf-core/nanoplot/main.nf'
include { HOSTILE_CLEAN } from '../../modules/nf-core/hostile/clean/main.nf'
include { CHOPPER } from '../../modules/nf-core/chopper/main.nf'


// -------------------------
// Workflow Definition
// -------------------------

workflow PREPROCESSING {

    take:
    reads_ch // channel: [ meta, reads ]
    ch_hostile_db // channel: tuple [db_name, db_dir]
    contam_ref // channel: path to contaminant reference

    main:

    ch_versions = channel.empty() // collect versions from all modules
    
    RawNanoPlot(reads_ch) // initial QC
    ch_versions = ch_versions.mix(RawNanoPlot.out.versions.first())
    
    // if using a host derived sample (human, mouse supported for now)
    if (params.run_hostile) {
        HOSTILE_CLEAN(reads_ch, ch_hostile_db) // host contamination removal
        ch_versions = ch_versions.mix(HOSTILE_CLEAN.out.versions)
        reads_dehosted_ch = HOSTILE_CLEAN.out.fastq
    } else {
        reads_dehosted_ch = reads_ch
    }

    CHOPPER(reads_dehosted_ch, contam_ref) // quality filtering; future contam removal
    ch_versions = ch_versions.mix(CHOPPER.out.versions.first())
    
    FinalNanoPlot(CHOPPER.out.fastq) // final QC

    emit:
    clean_reads = CHOPPER.out.fastq
    versions = ch_versions
}
