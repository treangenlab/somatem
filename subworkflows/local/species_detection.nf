#!/usr/bin/env nextflow

include { SYLPH_PROFILE } from "../../modules/nf-core/sylph/profile/main.nf"
include { GANON_CLASSIFY } from "../../modules/nf-core/ganon/classify/main.nf"
include { KRAKEN2_KRAKEN2 } from "../../modules/nf-core/kraken2/kraken2/main.nf"

// -------------------------
// Parameters
// -------------------------


// -------------------------
// Workflow
// -------------------------
workflow SPECIES_DETECTION {
    
    take:
    clean_reads_ch
    
    main: 
    
    taxonomy_report = channel.empty() // collect taxonomy reports from all tools

    SYLPH_PROFILE(clean_reads_ch, params.sylph_db)
    taxonomy_report = taxonomy_report.mix(SYLPH_PROFILE.out.profile_out)

    GANON_CLASSIFY(clean_reads_ch, params.ganon_db)
    taxonomy_report = taxonomy_report.mix(GANON_CLASSIFY.out.report)

    KRAKEN2_KRAKEN2(clean_reads_ch, params.kraken2_db, 'true', 'true')
    taxonomy_report = taxonomy_report.mix(KRAKEN2_KRAKEN2.out.report)
    
    emit:
    taxonomy_report = taxonomy_report
}