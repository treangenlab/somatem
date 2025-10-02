/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_somatem_pipeline'
include { DOWNLOAD_DBS } from '../subworkflows/local/download_dbs.nf'
include { PREPROCESSING } from '../subworkflows/local/pre-processing.nf'
include { TAXONOMIC_PROFILING } from '../subworkflows/local/taxonomic-profiling.nf'
include { GENOME_DYNAMICS } from '../subworkflows/local/genome-dynamics.nf'
include { ASSEMBLY_MAGS } from '../subworkflows/local/assembly_mags.nf'

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
    DOWNLOAD_DBS(0)


    // -----------------------------------------------------------------
    // Pre-processing and quality control on raw reads
    // -----------------------------------------------------------------
    contam_ref = Channel.value([]) // empty channel for now
    PREPROCESSING(ch_samplesheet, contam_ref)
    ch_versions = ch_versions.mix(PREPROCESSING.out.versions)

    // -----------------------------------------------------------------
    // Taxonomic profiling
    // -----------------------------------------------------------------
    
    if (params.analysis_type == "taxonomic-profiling") {
        TAXONOMIC_PROFILING(PREPROCESSING.out.clean_reads)
        ch_versions = ch_versions.mix(TAXONOMIC_PROFILING.out.versions)
        
        ch_key_outputs = ch_key_outputs.mix(TAXONOMIC_PROFILING.out.taxonomy_report)
    }

    // -----------------------------------------------------------------
    // assembly
    // -----------------------------------------------------------------
    if (params.analysis_type == "assembly") {

        // unpack the downloaded databases
        ch_checkm2_db = DOWNLOAD_DBS.out.ch_checkm2_db.map { _meta, db -> db } // strip meta, only take db
        ch_bakta_db = DOWNLOAD_DBS.out.ch_bakta_db
        ch_singlem_db = DOWNLOAD_DBS.out.ch_singlem_db

        ASSEMBLY_MAGS(PREPROCESSING.out.clean_reads, 
                ch_checkm2_db,
                ch_bakta_db,
                ch_singlem_db
        )
    }

    // -----------------------------------------------------------------
    // genome dynamics : Longitudinal analysis
    // -----------------------------------------------------------------
    if (params.analysis_type == "genome-dynamics") {
        GENOME_DYNAMICS(PREPROCESSING.out.clean_reads)
        ch_versions = ch_versions.mix(GENOME_DYNAMICS.out.versions)
        ch_key_outputs = ch_key_outputs.mix(GENOME_DYNAMICS.out.assembly_graph)
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
        ).set { ch_collated_versions }

    emit:
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
    clean_reads    = PREPROCESSING.out.clean_reads
    key_outputs    = ch_key_outputs              // channel: [ path(taxonomy_report.tsv) | path(assembly_graph.gfa), path(bandage_image.png) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
