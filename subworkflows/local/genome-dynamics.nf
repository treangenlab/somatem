#!/usr/bin/env nextflow

include { RHEA } from '../../modules/local/rhea/main.nf'
include { BANDAGE_IMAGE } from '../../modules/nf-core/bandage/image/main.nf'

// -------------------------
// Workflow Definition
// -------------------------
workflow GENOME_DYNAMICS {

    take:
    clean_reads_ch

    main:

    ch_versions = Channel.empty() // collect versions from all modules

    // collect reads and discard meta from clean_reads_ch
    collected_reads_ch = clean_reads_ch.map { meta, reads -> reads }
        .collect() // collect all reads into a single channel
        .map { reads -> 
            def flattened_reads = reads.flatten()
            return [[id:"multiple", single_end:true], flattened_reads] // add a mock meta
        }

    // debug
    collected_reads_ch.view()   

    // run Rhea
    RHEA(collected_reads_ch)

    // collect versions
    ch_versions = ch_versions.mix(RHEA.out.versions)

    // visualise assembly graph
    BANDAGE_IMAGE(RHEA.out.assembly_graph)
    ch_versions = ch_versions.mix(BANDAGE_IMAGE.out.versions)

    // collect outputs
    assembly_graph_outputs = RHEA.out.assembly_graph.mix(BANDAGE_IMAGE.out.png)

    emit:
    assembly_graph     = assembly_graph_outputs             // channel: [ path(assembly_graph.gfa), path(bandage_image.png) ]
    versions           = ch_versions                        // channel: [ path(versions.yml) ]
}