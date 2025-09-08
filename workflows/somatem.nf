/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_somatem_pipeline'
include { PREPROCESSING } from '../subworkflows/local/pre-processing.nf'
include { TAXONOMIC_PROFILING } from '../subworkflows/local/taxonomic-profiling.nf'

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

    // -----------------------------------------------------------------
    // Pre-processing and quality control on raw reads
    // -----------------------------------------------------------------
    contam_ref = Channel.value([]) // empty channel for now
    PREPROCESSING(ch_samplesheet, contam_ref)

    // -----------------------------------------------------------------
    // Taxonomic profiling
    // -----------------------------------------------------------------
    TAXONOMIC_PROFILING(PREPROCESSING.out.clean_reads)

    emit:
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
    clean_reads    = PREPROCESSING.out.clean_reads
    taxonomy_report = TAXONOMIC_PROFILING.out.taxonomy_report
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
