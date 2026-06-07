#!/usr/bin/env nextflow

include { RHEA } from '../../modules/local/rhea/main.nf'
include { BANDAGE_IMAGE } from '../../modules/nf-core/bandage/image/main.nf'
include { SOMATEM_SUMMARY_REPORT as GENOME_DYNAMICS_SUMMARY_REPORT } from '../../modules/local/somatem_summary_report/main.nf'

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

    // collect outputs
    assembly_graph_outputs = RHEA.out.assembly_graph.mix(BANDAGE_IMAGE.out.png)

    ch_versions_for_report = ch_versions
    ch_genome_dynamics_report_files = RHEA.out.assembly_graph
        .mix(
            RHEA.out.bandage_metadata,
            RHEA.out.bp_counts,
            RHEA.out.edge_coverage,
            RHEA.out.node_coverage_norm,
            RHEA.out.node_coverage,
            RHEA.out.rhea_log,
            RHEA.out.structural_variants,
            RHEA.out.gaf_files,
            BANDAGE_IMAGE.out.png,
            BANDAGE_IMAGE.out.svg,
            ch_versions_for_report
        )

    GENOME_DYNAMICS_SUMMARY_REPORT(
        'genome_dynamics',
        'Genome dynamics',
        Channel.fromPath(params.input),
        ch_genome_dynamics_report_files.flatMap { item ->
            def report_file = item
            if (item instanceof Collection && item.size() >= 2) {
                report_file = item[1]
            }
            if (report_file instanceof Collection) {
                return report_file
            }
            return [report_file]
        }.collect()
    )
    ch_versions = ch_versions.mix(GENOME_DYNAMICS_SUMMARY_REPORT.out.versions)

    emit:
    assembly_graph     = assembly_graph_outputs             // channel: [ path(assembly_graph.gfa), path(bandage_image.png) ]
    summary_report     = GENOME_DYNAMICS_SUMMARY_REPORT.out.html
    versions           = ch_versions                        // channel: [ path(versions.yml) ]
}
