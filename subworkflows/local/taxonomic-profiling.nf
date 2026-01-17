#!/usr/bin/env nextflow

include { EMU_ABUNDANCE } from '../../modules/local/emu/main.nf'
include { LEMUR } from '../../modules/local/lemur/main.nf'
include { MAGNET } from '../../modules/local/magnet/main.nf'
include { TAXBURST_CONVERT } from '../../modules/local/taxburst_convert/main.nf'
include { TAXBURST } from '../../modules/local/taxburst/main.nf'

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

    ch_lemur_db // channel: [ path(lemur_db) ]
    ch_sylph_db // channel: [ path(sylph_db) ]
    // ch_emu_db // channel: [ path(emu_db) ] // using the storeDir location supplied by ext.args for now (following TRANA/gms_16S format) 

    main:

    ch_versions = channel.empty() // collect versions from all modules

    // 16S amplicon reads
    if (params.data_type == "16S") {
        // single marker gene (16S) based tax profiling
        EMU_ABUNDANCE(clean_reads_ch)

        taxonomy_report = EMU_ABUNDANCE.out.report
        ch_versions = ch_versions.mix(EMU_ABUNDANCE.out.versions)

        // visualise taxonomy with Taxburst
        TAXBURST_CONVERT(taxonomy_report, 'emu') // Convert EMU output to Taxburst input
        TAXBURST(TAXBURST_CONVERT.out.converted, 'krona')
        ch_versions = ch_versions.mix(TAXBURST.out.versions)

    // metagenomic reads (default)
    } else {
        // multi-marker gene (16S + 18S + ITS) based tax profiling
        LEMUR(clean_reads_ch, ch_lemur_db)

        taxonomy_report = LEMUR.out.report
        classification_report = taxonomy_report
            .map { _meta, classification -> classification } // drop meta
        ch_versions = ch_versions.mix(LEMUR.out.versions)

        // visualise taxonomy with Taxburst
        TAXBURST_CONVERT(taxonomy_report, 'lemur') // Convert LEMUR output to Taxburst input
        TAXBURST(TAXBURST_CONVERT.out.converted, 'krona')
        ch_versions = ch_versions.mix(TAXBURST.out.versions)

        // Correct false positives for low abundance taxa / low coverage
        MAGNET(clean_reads_ch, classification_report)
        
        taxonomy_report = taxonomy_report.mix(MAGNET.out.report)
        ch_versions = ch_versions.mix(MAGNET.out.versions)
        // TODO: make magnet conditional on `validate_presence_absence` or `polished_profile` param

    }

    emit:
    taxonomy_report
    versions           = ch_versions                        // channel: [ path(versions.yml) ]
}
