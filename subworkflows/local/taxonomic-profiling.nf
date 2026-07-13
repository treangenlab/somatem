#!/usr/bin/env nextflow

include { EMU_ABUNDANCE } from '../../modules/local/emu/main.nf'
include { LEMUR } from '../../modules/local/lemur/main.nf'
include { MAGNET } from '../../modules/local/magnet/main.nf'
include { SYLPH_PROFILE } from '../../modules/nf-core/sylph/profile/main.nf'
include { KRAKEN2_KRAKEN2 } from '../../modules/nf-core/kraken2/kraken2/main.nf'
include { SINGLEM_PIPE } from '../../modules/local/singlem/pipe/main.nf'
include { TAXBURST_CONVERT } from '../../modules/local/taxburst_convert/main.nf'
include { TAXBURST } from '../../modules/local/taxburst/main.nf'
include { SOMATEM_SUMMARY_REPORT as TAXONOMY_SUMMARY_REPORT } from '../../modules/local/somatem_summary_report/main.nf'

// -------------------------
// Parameters
// -------------------------
params.input_dir   = 'examples/data'


// -------------------------
// Workflow Definition
// -------------------------

workflow TAXONOMIC_PROFILING {

    take:
    clean_reads_ch

    ch_taxonomy_db // channel: database selected for the active profiler
    // ch_emu_db // channel: [ path(emu_db) ] // using the storeDir location supplied by ext.args for now (following TRANA/gms_16S format) 

    main:

    ch_versions = channel.empty() // collect versions from all modules
    taxonomy_visualization = channel.empty()
    taxonomy_report_slug = params.data_type == "16S" ? "16S" : "taxonomic_profiling"
    taxonomy_report_label = params.data_type == "16S" ? "16S taxonomic profiling" : "Metagenomic taxonomic profiling"

    // 16S amplicon reads
    if (params.data_type == "16S") {
        // single marker gene (16S) based tax profiling
        EMU_ABUNDANCE(clean_reads_ch)

        taxonomy_report = EMU_ABUNDANCE.out.report
        ch_versions = ch_versions.mix(EMU_ABUNDANCE.out.versions)

        // visualise taxonomy with Taxburst
        TAXBURST_CONVERT(taxonomy_report, 'emu') // Convert EMU output to Taxburst input
        TAXBURST(TAXBURST_CONVERT.out.converted, 'krona')
        taxonomy_visualization = TAXBURST.out.html
        ch_versions = ch_versions.mix(TAXBURST.out.versions)

    // metagenomic reads (default)
    } else {
        if (params.taxonomic_profiler == 'sylph') {
            SYLPH_PROFILE(clean_reads_ch, ch_taxonomy_db)
            taxonomy_report = SYLPH_PROFILE.out.profile_out
            ch_versions = ch_versions.mix(SYLPH_PROFILE.out.versions)

            TAXBURST_CONVERT(taxonomy_report, 'sylph')
            TAXBURST(TAXBURST_CONVERT.out.converted, 'krona')
            taxonomy_visualization = TAXBURST.out.html
            ch_versions = ch_versions.mix(TAXBURST_CONVERT.out.versions, TAXBURST.out.versions)
        } else if (params.taxonomic_profiler == 'lemur-magnet') {
            LEMUR(clean_reads_ch, ch_taxonomy_db)
            taxonomy_report = LEMUR.out.report
            classification_report = taxonomy_report.map { _meta, classification -> classification }
            ch_versions = ch_versions.mix(LEMUR.out.versions)

            TAXBURST_CONVERT(taxonomy_report, 'lemur')
            TAXBURST(TAXBURST_CONVERT.out.converted, 'krona')
            taxonomy_visualization = TAXBURST.out.html
            ch_versions = ch_versions.mix(TAXBURST.out.versions)

            MAGNET(clean_reads_ch, classification_report)
            taxonomy_report = taxonomy_report.mix(MAGNET.out.report)
            ch_versions = ch_versions.mix(MAGNET.out.versions)
        } else if (params.taxonomic_profiler == 'kraken2') {
            KRAKEN2_KRAKEN2(clean_reads_ch, ch_taxonomy_db, params.kraken2_save_output_fastqs, params.kraken2_save_reads_assignment)
            taxonomy_report = KRAKEN2_KRAKEN2.out.report
            ch_versions = ch_versions.mix(KRAKEN2_KRAKEN2.out.versions)
        } else if (params.taxonomic_profiler == 'singlem') {
            SINGLEM_PIPE(clean_reads_ch, ch_taxonomy_db, 'reads')
            taxonomy_report = SINGLEM_PIPE.out.taxonomic_profile
            taxonomy_visualization = SINGLEM_PIPE.out.krona_profile
            ch_versions = ch_versions.mix(SINGLEM_PIPE.out.versions)
        } else {
            error("Unsupported taxonomic profiler '${params.taxonomic_profiler}'. Choose sylph, lemur-magnet, kraken2, or singlem.")
        }

    }

    ch_versions_for_report = ch_versions
    ch_taxonomy_report_files = taxonomy_report.mix(taxonomy_visualization, ch_versions_for_report)
    TAXONOMY_SUMMARY_REPORT(
        taxonomy_report_slug,
        taxonomy_report_label,
        channel.fromPath(params.input),
        ch_taxonomy_report_files.flatMap { item ->
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
    ch_versions = ch_versions.mix(TAXONOMY_SUMMARY_REPORT.out.versions)

    emit:
    taxonomy_report
    summary_report      = TAXONOMY_SUMMARY_REPORT.out.html
    versions           = ch_versions                        // channel: [ path(versions.yml) ]
}
