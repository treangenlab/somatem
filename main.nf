#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/somatem
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/somatem
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SOMATEM  } from './workflows/somatem'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_somatem_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_somatem_pipeline'
include { getGenomeAttribute      } from './subworkflows/local/utils_nfcore_somatem_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// TODO nf-core: Remove this line if you don't need a FASTA file
//   This is an example of how to use getGenomeAttribute() to fetch parameters
//   from igenomes.config using `--genome`
params.fasta = getGenomeAttribute('fasta')

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow ORCHESTRATE_SOMATEM {

    take:
    samplesheet // channel: samplesheet read in from --input

    main:

    //
    // WORKFLOW: Run pipeline
    //
    SOMATEM (
        samplesheet
    )

    emit:
    versions        = SOMATEM.out.versions                 // channel: [ path(versions.yml) ]
    clean_reads     = SOMATEM.out.clean_reads
    key_outputs     = SOMATEM.out.key_outputs              // channel: [ path(taxonomy_report.tsv) | path(assembly_graph.gfa), path(bandage_image.png) ]

    // Assembly outputs
    mapping         = SOMATEM.out.mapping                  // channel: [ val(meta), path(*.bam) ]
    bin_tables      = SOMATEM.out.bin_tables               // channel: [ path(*.csv) | path(*.tsv) ]
    bin_fasta       = SOMATEM.out.bin_fasta                // channel: [ path(*.fa.gz) ]
    assembly        = SOMATEM.out.assembly                 // channel: [ val(meta), path(*.fasta.gz) ]
    assembly_gfa    = SOMATEM.out.assembly_gfa             // channel: [ val(meta), path(*.gfa.gz) ]
    assembly_log    = SOMATEM.out.assembly_log             // channel: [ val(meta), path(*.log) ]
    coverage        = SOMATEM.out.coverage                 // channel: [ val(meta), path(*.txt) ]
    checkm2_report  = SOMATEM.out.checkm2_report           // channel: [ val(meta), path(*.tsv) ]
    
    // Annotation outputs
    bakta_gff       = SOMATEM.out.bakta_gff                // channel: [ val(meta), path(*.gff) ]
    bakta_tsv       = SOMATEM.out.bakta_tsv                // channel: [ val(meta), path(*.tsv) ]
    bakta_txt       = SOMATEM.out.bakta_txt                // channel: [ val(meta), path(*.txt) ]
    
    // Taxonomic outputs
    singlem_profile = SOMATEM.out.singlem_profile          // channel: [ val(meta), path(*.tsv) ]
    taxburst_html   = SOMATEM.out.taxburst_html            // channel: [ val(meta), path(*.html) ]
    
    // Analysis outputs
    pigeon_html     = SOMATEM.out.pigeon_html              // channel: [ val(meta), path(*.html) ]
    appraise_summary = SOMATEM.out.appraise_summary        // channel: [ val(meta), path(*.tsv) ]

}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GLOBAL OUTPUT DEFINITIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

output {
    // Quality Control outputs
    clean_reads {
        path { meta, _reads -> "preprocessing/${meta.id}" }
    }
    
    // Assembly outputs
    assembly_fasta {
        path { meta, _fasta -> "assembly/${meta.id}" }
    }
    assembly_gfa {
        path { meta, _gfa -> "assembly/${meta.id}" }
    }
    assembly_logs {
        path { meta, _log -> "assembly/${meta.id}/logs" }
    }
    
    // Mapping outputs
    mapping {
        path { meta, _bam -> "mapping/${meta.id}" }
    }
    coverage {
        path { meta, _cov -> "mapping/${meta.id}/coverage" }
    }
    
    // Binning outputs
    binning_tables {
        path { "binning/tables" }
    }
    binning_fasta {
        path { "binning/fasta" }
    }
    
    // Quality assessment
    checkm2_reports {
        path { meta, _report -> "quality_assessment/${meta.id}" }
    }
    
    // Annotation outputs (high-quality bins)
    annotation_gff {
        path { meta, _gff -> "annotation/${meta.sample_id}/${meta.id}" }
    }
    annotation_tables {
        path { meta, _tsv -> "annotation/${meta.sample_id}/${meta.id}" }
    }
    annotation_summary {
        path { meta, _txt -> "annotation/${meta.sample_id}/${meta.id}" }
    }
    
    // Taxonomic profiling outputs
    taxonomy_profiles {
        path { meta, _profile -> "taxonomy/${meta.id}" }
    }
    taxonomy_visualizations {
        path { meta, _html -> "taxonomy/${meta.id}/visualizations" }
    }
    
    // Post-hoc analysis
    pigeon_reports {
        path { meta, _html -> "post_hoc_analysis/${meta.id}/pigeon" }
    }
    appraise_reports {
        path { meta, _summary -> "post_hoc_analysis/${meta.id}/appraise" }
    }
    
    // Pipeline metadata
    versions {
        path { "pipeline_info" }
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input
    )

    //
    // WORKFLOW: Run main workflow
    //
    ORCHESTRATE_SOMATEM (
        PIPELINE_INITIALISATION.out.samplesheet
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.outdir,
        params.monochrome_logs,
    )

    publish:
    // Quality control outputs
    clean_reads = ORCHESTRATE_SOMATEM.out.clean_reads
    
    // Assembly outputs
    assembly_fasta = ORCHESTRATE_SOMATEM.out.assembly
    assembly_gfa = ORCHESTRATE_SOMATEM.out.assembly_gfa
    assembly_logs = ORCHESTRATE_SOMATEM.out.assembly_log
    
    // Mapping outputs
    mapping = ORCHESTRATE_SOMATEM.out.mapping
    coverage = ORCHESTRATE_SOMATEM.out.coverage
    
    // Binning outputs  
    binning_tables = ORCHESTRATE_SOMATEM.out.bin_tables
    binning_fasta = ORCHESTRATE_SOMATEM.out.bin_fasta
    
    // Quality assessment
    checkm2_reports = ORCHESTRATE_SOMATEM.out.checkm2_report
    
    // Annotation outputs (high-quality bins)
    annotation_gff = ORCHESTRATE_SOMATEM.out.bakta_gff
    annotation_tables = ORCHESTRATE_SOMATEM.out.bakta_tsv
    annotation_summary = ORCHESTRATE_SOMATEM.out.bakta_txt
    
    // Taxonomic profiling outputs
    taxonomy_profiles = ORCHESTRATE_SOMATEM.out.singlem_profile
    taxonomy_visualizations = ORCHESTRATE_SOMATEM.out.taxburst_html
    
    // Post-hoc analysis
    pigeon_reports = ORCHESTRATE_SOMATEM.out.pigeon_html
    appraise_reports = ORCHESTRATE_SOMATEM.out.appraise_summary
    
    // Version tracking
    versions = ORCHESTRATE_SOMATEM.out.versions
    
    onComplete:
    log.info "=" * 80
    log.info "Pipeline execution completed!"
    log.info "Status: ${workflow.success ? 'SUCCESS' : 'FAILED'}"
    log.info "Completed at: ${workflow.complete}"
    log.info "Duration: ${workflow.duration}"
    log.info "Results published to: ${params.outdir}"
    log.info "=" * 80
    
    onError:
    log.error "=" * 80
    log.error "Pipeline execution failed!"
    log.error "Error message: ${workflow.errorMessage}"
    log.error "Failed at: ${workflow.errorReport}"
    log.error "=" * 80
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
