/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_somatem_pipeline'
include { DOWNLOAD_DBS } from '../subworkflows/local/download_databases.nf'
include { PREPROCESSING } from '../subworkflows/local/pre-processing.nf'
include { TAXONOMIC_PROFILING } from '../subworkflows/local/taxonomic-profiling.nf'
include { GENOME_DYNAMICS } from '../subworkflows/local/genome-dynamics.nf'
include { ASSEMBLY_MAGS } from '../subworkflows/local/assembly_mags.nf'
include { ISOLATE_ANALYSIS } from '../subworkflows/local/isolate_analysis.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SOMATEM {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    ch_versions = Channel.empty()
    ch_key_outputs = Channel.empty()

    // -----------------------------------------------------------------
    // Download databases
    // -----------------------------------------------------------------
    DOWNLOAD_DBS(params.analysis_type, params.hostile_index, 
            params.lemur_db_zenodo_id, params.checkm2_db_zenodo_id)


    // -----------------------------------------------------------------
    // Pre-processing and quality control on raw reads
    // -----------------------------------------------------------------
    contam_ref = Channel.value([]) // empty channel for now
    if (params.analysis_type == "isolate-analysis") {
        ch_clean_reads = ch_samplesheet
        ch_summary_reports = Channel.empty()
    } else {
        PREPROCESSING(ch_samplesheet, DOWNLOAD_DBS.out.ch_hostile_db, contam_ref)
        ch_versions = ch_versions.mix(PREPROCESSING.out.versions)
        ch_clean_reads = PREPROCESSING.out.clean_reads
        ch_summary_reports = PREPROCESSING.out.summary_report
    }

    // -----------------------------------------------------------------
    // Taxonomic profiling
    // -----------------------------------------------------------------
    
    if (params.analysis_type == "taxonomic-profiling") {
        TAXONOMIC_PROFILING(ch_clean_reads, DOWNLOAD_DBS.out.ch_lemur_db)
        ch_versions = ch_versions.mix(TAXONOMIC_PROFILING.out.versions)
        
        ch_key_outputs = ch_key_outputs.mix(TAXONOMIC_PROFILING.out.taxonomy_report)
        ch_summary_reports = ch_summary_reports.mix(TAXONOMIC_PROFILING.out.summary_report)
    }

    // -----------------------------------------------------------------
    // assembly
    // -----------------------------------------------------------------
    if (params.analysis_type == "assembly") {

        // unpack the downloaded databases
        ch_checkm2_db = DOWNLOAD_DBS.out.ch_checkm2_db.map { _meta, db -> db } // strip meta, only take db
        ch_bakta_db = DOWNLOAD_DBS.out.ch_bakta_db
        ch_singlem_db = DOWNLOAD_DBS.out.ch_singlem_db

        ASSEMBLY_MAGS(ch_clean_reads, 
                ch_checkm2_db,
                ch_bakta_db,
                ch_singlem_db
        )
        ch_versions = ch_versions.mix(ASSEMBLY_MAGS.out.versions)
        ch_summary_reports = ch_summary_reports.mix(ASSEMBLY_MAGS.out.summary_report)

        // collect key outputs: not using right now ; have separate emits below

    }

    // -----------------------------------------------------------------
    // isolate analysis
    // -----------------------------------------------------------------
    if (params.analysis_type == "isolate-analysis") {

        ch_bakta_db = DOWNLOAD_DBS.out.ch_bakta_db
        ch_kraken2_db = DOWNLOAD_DBS.out.ch_kraken2_db
        ch_checkm2_db = DOWNLOAD_DBS.out.ch_checkm2_db.map { _meta, db -> db }

        ISOLATE_ANALYSIS(
            ch_clean_reads,
            ch_bakta_db,
            ch_kraken2_db,
            ch_checkm2_db
        )
        ch_versions = ch_versions.mix(ISOLATE_ANALYSIS.out.versions)
        ch_key_outputs = ch_key_outputs.mix(
            ISOLATE_ANALYSIS.out.assembly,
            ISOLATE_ANALYSIS.out.kraken2_report,
            ISOLATE_ANALYSIS.out.bakta_tsv,
            ISOLATE_ANALYSIS.out.checkm2_tsv
        )
        ch_key_outputs = ch_key_outputs.mix(ISOLATE_ANALYSIS.out.btyper3_results)
        ch_summary_reports = ch_summary_reports.mix(ISOLATE_ANALYSIS.out.summary_report)
    }

    // -----------------------------------------------------------------
    // genome dynamics : Longitudinal analysis
    // -----------------------------------------------------------------
    if (params.analysis_type == "genome-dynamics") {
        GENOME_DYNAMICS(ch_clean_reads)
        ch_versions = ch_versions.mix(GENOME_DYNAMICS.out.versions)
        ch_key_outputs = ch_key_outputs.mix(GENOME_DYNAMICS.out.assembly_graph)
        ch_summary_reports = ch_summary_reports.mix(GENOME_DYNAMICS.out.summary_report)
    }


    // -----------------------------------------------------------------
    // Collate and save software versions
    // -----------------------------------------------------------------
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'somatem_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { _ch_collated_versions }

    emit:
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
    clean_reads    = ch_clean_reads
    summary_reports = ch_summary_reports
    key_outputs    = ch_key_outputs              // channel: [ path(taxonomy_report.tsv) | path(assembly_graph.gfa), path(bandage_image.png) ]
    
    // separate key emits for publishing convenience
    mapping       = params.analysis_type == "assembly" ? ASSEMBLY_MAGS.out.bam_sorted : channel.empty()  // channel: [ val(meta), path(*.bam) ]
    bin_tables    = params.analysis_type == "assembly" ? ASSEMBLY_MAGS.out.bins_csv.mix(ASSEMBLY_MAGS.out.bins_tsv, ASSEMBLY_MAGS.out.iterative_manifest, ASSEMBLY_MAGS.out.iterative_selected, ASSEMBLY_MAGS.out.iterative_trajectory, ASSEMBLY_MAGS.out.iterative_summary, ASSEMBLY_MAGS.out.iterative_command_log) : channel.empty() // channel: [ path(*.csv) | path(*.tsv) ]
    bin_fasta     = params.analysis_type == "assembly" ? ASSEMBLY_MAGS.out.bins : channel.empty() // channel: [ path(*.fa.gz) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
