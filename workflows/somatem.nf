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
    PREPROCESSING(ch_samplesheet, DOWNLOAD_DBS.out.ch_hostile_db, contam_ref)
    ch_versions = ch_versions.mix(PREPROCESSING.out.versions)

    // -----------------------------------------------------------------
    // Taxonomic profiling
    // -----------------------------------------------------------------
    
    if (params.analysis_type == "taxonomic-profiling") {
        TAXONOMIC_PROFILING(PREPROCESSING.out.clean_reads, DOWNLOAD_DBS.out.ch_lemur_db)
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
        ch_versions = ch_versions.mix(ASSEMBLY_MAGS.out.versions)

        // collect key outputs: not using right now ; have separate emits below

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
        ).set { _ch_collated_versions }

    emit:
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
    clean_reads    = PREPROCESSING.out.clean_reads
    key_outputs    = ch_key_outputs              // channel: [ path(taxonomy_report.tsv) | path(assembly_graph.gfa), path(bandage_image.png) ]
    
    // Assembly outputs
    mapping        = params.analysis_type == "assembly" ? ASSEMBLY_MAGS.out.bam_sorted : channel.empty()  // channel: [ val(meta), path(*.bam) ]
    bin_tables     = params.analysis_type == "assembly" ? ASSEMBLY_MAGS.out.bins_csv.mix(ASSEMBLY_MAGS.out.bins_tsv) : channel.empty() // channel: [ path(*.csv) | path(*.tsv) ]
    bin_fasta      = params.analysis_type == "assembly" ? ASSEMBLY_MAGS.out.bins : channel.empty() // channel: [ path(*.fa.gz) ]
    assembly       = params.analysis_type == "assembly" ? ASSEMBLY_MAGS.out.assembly : channel.empty() // channel: [ val(meta), path(*.fasta.gz) ]
    assembly_gfa   = params.analysis_type == "assembly" ? ASSEMBLY_MAGS.out.assembly_gfa : channel.empty() // channel: [ val(meta), path(*.gfa.gz) ]
    assembly_log   = params.analysis_type == "assembly" ? ASSEMBLY_MAGS.out.assembly_log : channel.empty() // channel: [ val(meta), path(*.log) ]
    coverage       = params.analysis_type == "assembly" ? ASSEMBLY_MAGS.out.coverage : channel.empty() // channel: [ val(meta), path(*.txt) ]
    checkm2_report = params.analysis_type == "assembly" ? ASSEMBLY_MAGS.out.checkm2_report : channel.empty() // channel: [ val(meta), path(*.tsv) ]
    
    // Annotation outputs (high-quality bins only)
    bakta_gff      = params.analysis_type == "assembly" ? ASSEMBLY_MAGS.out.bakta_gff : channel.empty() // channel: [ val(meta), path(*.gff) ]
    bakta_tsv      = params.analysis_type == "assembly" ? ASSEMBLY_MAGS.out.bakta_tsv : channel.empty() // channel: [ val(meta), path(*.tsv) ]
    bakta_txt      = params.analysis_type == "assembly" ? ASSEMBLY_MAGS.out.bakta_txt : channel.empty() // channel: [ val(meta), path(*.txt) ]
    
    // Taxonomic profiling outputs
    singlem_profile = params.analysis_type == "assembly" ? ASSEMBLY_MAGS.out.singlem_profile : channel.empty() // channel: [ val(meta), path(*.tsv) ]
    taxburst_html   = params.analysis_type == "assembly" ? ASSEMBLY_MAGS.out.taxburst_html : channel.empty() // channel: [ val(meta), path(*.html) ]
    
    // Post-hoc analysis
    pigeon_html     = params.analysis_type == "assembly" ? ASSEMBLY_MAGS.out.pigeon_html : channel.empty() // channel: [ val(meta), path(*.html) ]
    appraise_summary = params.analysis_type == "assembly" ? ASSEMBLY_MAGS.out.appraise_summary : channel.empty() // channel: [ val(meta), path(*.tsv) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
