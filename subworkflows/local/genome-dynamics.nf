#!/usr/bin/env nextflow

include { RHEA } from '../../modules/local/rhea/main.nf'

// -------------------------
// Workflow Definition
// -------------------------
workflow GENOME_DYNAMICS {

    take:
    clean_reads_ch

    main:

    ch_versions = Channel.empty() // collect versions from all modules

    // run Rhea
    RHEA(clean_reads_ch)

    // rhea_out = RHEA.out.some_channel
    ch_versions = ch_versions.mix(RHEA.out.versions)

    // visualise assembly graph
        

    emit:
    versions           = ch_versions                        // channel: [ path(versions.yml) ]
}