#!/usr/bin/env nextflow

include { convert_to_nfcore_tuple } from './utils/nf-core-compatibility.nf'

include { NANOPLOT as RawNanoPlot; NANOPLOT as FinalNanoPlot } from '../../modules/nf-core/nanoplot/main.nf'
include { HOSTILE_CLEAN } from '../../modules/nf-core/hostile/clean/main.nf'
include { CHOPPER } from '../../modules/nf-core/chopper/main.nf'
include { SOMATEM_SUMMARY_REPORT as PREPROCESSING_SUMMARY_REPORT } from '../../modules/local/somatem_summary_report/main.nf'


// -------------------------
// Workflow Definition
// -------------------------

workflow PREPROCESSING {

    take:
    reads_ch // channel: [ meta, reads ]
    ch_hostile_db // channel: tuple [db_name, db_dir]
    contam_ref // channel: path to contaminant reference

    main:

    ch_versions = channel.empty() // collect versions from all modules
    
    RawNanoPlot(reads_ch) // initial QC
    ch_versions = ch_versions.mix(RawNanoPlot.out.versions.first())
    
    // if using a host derived sample (human, mouse supported for now)
    if (params.run_hostile) {
        HOSTILE_CLEAN(reads_ch, ch_hostile_db) // host contamination removal
        ch_versions = ch_versions.mix(HOSTILE_CLEAN.out.versions)
        reads_dehosted_ch = HOSTILE_CLEAN.out.fastq
    } else {
        reads_dehosted_ch = reads_ch
    }

    CHOPPER(reads_dehosted_ch, contam_ref) // quality filtering; future contam removal
    ch_versions = ch_versions.mix(CHOPPER.out.versions.first())
    
    FinalNanoPlot(CHOPPER.out.fastq) // final QC
    ch_versions = ch_versions.mix(FinalNanoPlot.out.versions.first())

    ch_versions_for_report = ch_versions
    ch_preprocessing_report_files = RawNanoPlot.out.txt
        .mix(FinalNanoPlot.out.txt, CHOPPER.out.fastq)
        .mix(ch_versions_for_report)

    if (params.run_hostile) {
        ch_preprocessing_report_files = ch_preprocessing_report_files.mix(HOSTILE_CLEAN.out.fastq)
    }

    PREPROCESSING_SUMMARY_REPORT(
        'pre_processing',
        'Pre-processing and quality control',
        Channel.fromPath(params.input),
        ch_preprocessing_report_files.flatMap { item ->
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
    ch_versions = ch_versions.mix(PREPROCESSING_SUMMARY_REPORT.out.versions)

    emit:
    clean_reads = CHOPPER.out.fastq
    summary_report = PREPROCESSING_SUMMARY_REPORT.out.html
    versions = ch_versions
}
