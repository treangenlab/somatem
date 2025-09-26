/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_somatem_pipeline'
include { PREPROCESSING } from '../subworkflows/local/pre-processing.nf'
include { TAXONOMIC_PROFILING } from '../subworkflows/local/taxonomic-profiling.nf'
include { GENOME_DYNAMICS } from '../subworkflows/local/genome-dynamics.nf'

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
    }

    // -----------------------------------------------------------------
    // assembly
    // -----------------------------------------------------------------
    // if (params.analysis_type == "assembly") {
    //     ASSEMBLY(clean_reads)
    // }

    // -----------------------------------------------------------------
    // genome dynamics : Longitudinal analysis
    // -----------------------------------------------------------------
    if (params.analysis_type == "genome-dynamics") {
        GENOME_DYNAMICS(PREPROCESSING.out.clean_reads)
        ch_versions = ch_versions.mix(GENOME_DYNAMICS.out.versions)
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
    // taxonomy_report = TAXONOMIC_PROFILING.out.taxonomy_report
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
