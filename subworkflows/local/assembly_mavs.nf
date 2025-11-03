#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PARAMETERS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params.input = null
params.outdir = './results'
params.host_fasta = null
params.min_len = 1500
params.run_vcontact2 = true

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// modules
include { FLYE                   } from '../../modules/nf-core/flye/main'
include { VIRSORTER2             } from '../../modules/local/virsorter2/main'
include { CHECKV                 } from '../../modules/local/checkv/main'
include { COVERM                 } from '../../modules/local/coverm/main'
include { VRHYME                 } from '../../modules/local/vrhyme/main'
include { IPHOP                  } from '../../modules/local/iphop/main'
include { VCONTACT2              } from '../../modules/local/vcontact2/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: VIRAL METAGENOMICS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ASSEMBLY_MAVS {

    take:
    ch_clean_reads    // channel: [ val(meta), path(reads) ]
    ch_host_fasta     // channel: [ path(host_fasta) ] - optional

    main:
    ch_versions = Channel.empty()

    //
    // Long-read assembly with metaFlye
    //
    FLYE(ch_clean_reads, params.flye_mode)
    ch_versions = ch_versions.mix(FLYE.out.versions.first())

    //
    // Viral identification with VirSorter2
    //
    VIRSORTER2(
        FLYE.out.fasta,
        params.min_len ?: 1500
    )
    ch_versions = ch_versions.mix(VIRSORTER2.out.versions.first())

    //
    // Quality control and provirus trimming with CheckV
    //
    CHECKV(VIRSORTER2.out.viral)
    ch_versions = ch_versions.mix(CHECKV.out.versions.first())

    //
    // Coverage analysis with CoverM
    //
    COVERM(
        ch_clean_reads,
        CHECKV.out.curated
    )
    ch_versions = ch_versions.mix(COVERM.out.versions.first())

    //
    // Viral binning with vRhyme
    //
    VRHYME(
        CHECKV.out.curated,
        COVERM.out.coverage
    )
    ch_versions = ch_versions.mix(VRHYME.out.versions.first())

    //
    // Host prediction with iPHoP
    //
    IPHOP(VRHYME.out.bins)
    ch_versions = ch_versions.mix(IPHOP.out.versions.first())

    //
    // Taxonomy with vConTACT2
    //
    ch_taxonomy_results = Channel.empty()
    if (params.run_vcontact2 ?: true) {
        VCONTACT2(VRHYME.out.bins)
        ch_versions = ch_versions.mix(VCONTACT2.out.versions.first())
        ch_taxonomy_results = VCONTACT2.out.overview
    }

    emit:
    versions            = ch_versions                         // channel: [ path(versions.yml) ]
    vmags              = VRHYME.out.bins                      // channel: [ val(meta), path(bins) ]
    curated_contigs    = CHECKV.out.curated                   // channel: [ val(meta), path(contigs) ]
    host_predictions   = IPHOP.out.genus_predictions          // channel: [ val(meta), path(predictions) ]
    taxonomy_results   = ch_taxonomy_results                  // channel: [ val(meta), path(taxonomy) ]
    assembly_stats     = FLYE.out.log                         // channel: [ val(meta), path(log) ]
    viral_quality      = CHECKV.out.quality                   // channel: [ val(meta), path(quality) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MAIN ENTRY WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    
    // Validate required parameters
    if (!params.input) {
        error "Please provide an input samplesheet with --input"
    }

    // Create input channel from samplesheet
    ch_input = Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def meta = [id: row.sample]
            def reads = []
            
            // Handle paired-end reads
            if (row.fastq_1 && row.fastq_2) {
                reads = [file(row.fastq_1, checkIfExists: true), file(row.fastq_2, checkIfExists: true)]
            }
            // Handle single-end reads
            else if (row.fastq_1) {
                reads = [file(row.fastq_1, checkIfExists: true)]
            }
            else {
                error "No FASTQ files specified for sample ${row.sample}"
            }
            
            return [meta, reads]
        }

    // Create host fasta channel if provided
    ch_host_fasta = params.host_fasta ? 
        Channel.fromPath(params.host_fasta, checkIfExists: true) : 
        Channel.empty()

    // Run the ASSEMBLY_MAVS subworkflow
    ASSEMBLY_MAVS(ch_input, ch_host_fasta)

    // Publish results
    ASSEMBLY_MAVS.out.vmags
        .subscribe { meta, files ->
            files.each { file ->
                file.copyTo("${params.outdir}/vmags/${meta.id}_${file.name}")
            }
        }

    ASSEMBLY_MAVS.out.curated_contigs
        .subscribe { meta, files ->
            files.each { file ->
                file.copyTo("${params.outdir}/curated_contigs/${meta.id}_${file.name}")
            }
        }

    ASSEMBLY_MAVS.out.host_predictions
        .subscribe { meta, files ->
            files.each { file ->
                file.copyTo("${params.outdir}/host_predictions/${meta.id}_${file.name}")
            }
        }

    ASSEMBLY_MAVS.out.taxonomy_results
        .subscribe { meta, files ->
            files.each { file ->
                file.copyTo("${params.outdir}/taxonomy_results/${meta.id}_${file.name}")
            }
        }

    ASSEMBLY_MAVS.out.assembly_stats
        .subscribe { meta, files ->
            files.each { file ->
                file.copyTo("${params.outdir}/assembly_stats/${meta.id}_${file.name}")
            }
        }

    ASSEMBLY_MAVS.out.viral_quality
        .subscribe { meta, files ->
            files.each { file ->
                file.copyTo("${params.outdir}/viral_quality/${meta.id}_${file.name}")
            }
        }

    // Save versions
    ASSEMBLY_MAVS.out.versions
        .collectFile(name: 'software_versions.yml', storeDir: "${params.outdir}/pipeline_info")
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/