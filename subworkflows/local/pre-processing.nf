#!/usr/bin/env nextflow

include { convert_to_nfcore_tuple } from './utils/nf-core-compatibility.nf'

include { NANOPLOT as RawNanoPlot; NANOPLOT as FinalNanoPlot } from '../../modules/nf-core/nanoplot/main.nf'
include { DEACON_FILTER } from '../../modules/local/deacon/filter/main.nf'
include { HOSTILE_CLEAN } from '../../modules/local/hostile/clean/main.nf'
include { CHOPPER } from '../../modules/nf-core/chopper/main.nf'
include { SOMATEM_SUMMARY_REPORT as PREPROCESSING_SUMMARY_REPORT } from '../../modules/local/somatem_summary_report/main.nf'


// -------------------------
// Workflow Definition
// -------------------------

workflow PREPROCESSING {

    take:
    reads_ch // channel: [ meta, reads ]
    ch_host_index // reusable host index for the selected filter
    contam_ref // channel: path to contaminant reference

    main:

    ch_versions = channel.empty() // collect versions from all modules
    
    RawNanoPlot(reads_ch) // initial QC
    ch_versions = ch_versions.mix(RawNanoPlot.out.versions.first())
    
    ch_host_report_files = channel.empty()
    reads_dehosted_ch = reads_ch
    if (params.hostile_reference || params.run_deacon) {
        if (params.hostile_reference) {
            host_filter = HOSTILE_CLEAN(reads_ch, ch_host_index)
        } else {
            host_filter = DEACON_FILTER(reads_ch, ch_host_index)
        }
        reads_dehosted_ch = host_filter.fastq
        ch_versions = ch_versions.mix(host_filter.versions)
        ch_host_report_files = host_filter.fastq.mix(host_filter.json)
    }

    CHOPPER(reads_dehosted_ch, contam_ref) // quality filtering; future contam removal
    ch_versions = ch_versions.mix(CHOPPER.out.versions.first())
    
    FinalNanoPlot(CHOPPER.out.fastq) // final QC
    ch_versions = ch_versions.mix(FinalNanoPlot.out.versions.first())

    ch_versions_for_report = ch_versions
    ch_preprocessing_report_files = RawNanoPlot.out.txt
        .mix(FinalNanoPlot.out.txt, CHOPPER.out.fastq, ch_host_report_files)
        .mix(ch_versions_for_report)

    PREPROCESSING_SUMMARY_REPORT(
        'pre_processing',
        'Pre-processing and quality control',
        channel.fromPath(params.input),
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
