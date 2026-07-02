#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SeqScreen DSL2 subworkflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    This is a modularized port of the upstream SeqScreen short-read/contig pipeline.
    The upstream helper scripts are vendored under assets/seqscreen and staged into
    each module that needs them.
----------------------------------------------------------------------------------------
*/

include { SEQSCREEN_VALIDATE_FASTA } from '../../modules/local/seqscreen/validate_fasta/main'
include { SEQSCREEN_WINDOW_FASTA } from '../../modules/local/seqscreen/window_fasta/main'
include { SEQSCREEN_BOWTIE2 as SEQSCREEN_BOWTIE2_BSAT } from '../../modules/local/seqscreen/bowtie2/main'
include { SEQSCREEN_BOWTIE2 as SEQSCREEN_BOWTIE2_VFDB } from '../../modules/local/seqscreen/bowtie2/main'
include { SEQSCREEN_RAPSEARCH2 as SEQSCREEN_RAPSEARCH2_BSAT } from '../../modules/local/seqscreen/rapsearch2/main'
include { SEQSCREEN_RAPSEARCH2 as SEQSCREEN_RAPSEARCH2_VFDB } from '../../modules/local/seqscreen/rapsearch2/main'
include { SEQSCREEN_THREATS_BY_BLACKLIST as SEQSCREEN_THREATS_BY_BLACKLIST_BSAT } from '../../modules/local/seqscreen/threats_by_blacklist/main'
include { SEQSCREEN_THREATS_BY_BLACKLIST as SEQSCREEN_THREATS_BY_BLACKLIST_VFDB } from '../../modules/local/seqscreen/threats_by_blacklist/main'
include { SEQSCREEN_HMMSCAN } from '../../modules/local/seqscreen/hmmscan/main'
include { SEQSCREEN_MUMMER } from '../../modules/local/seqscreen/mummer/main'
include { SEQSCREEN_DIAMOND_BLASTX } from '../../modules/local/seqscreen/diamond_blastx/main'
include { SEQSCREEN_BLASTX } from '../../modules/local/seqscreen/blastx/main'
include { SEQSCREEN_CENTRIFUGE } from '../../modules/local/seqscreen/centrifuge/main'
include { SEQSCREEN_BLASTN as SEQSCREEN_BLASTN_NT } from '../../modules/local/seqscreen/blastn/main'
include { SEQSCREEN_BLASTN as SEQSCREEN_BLASTN_MEGARES } from '../../modules/local/seqscreen/blastn/main'
include { SEQSCREEN_OUTLIER_DETECTION } from '../../modules/local/seqscreen/outlier_detection/main'
include { SEQSCREEN_TAXONOMIC_ASSIGNMENT_FAST; SEQSCREEN_TAXONOMIC_ASSIGNMENT_SENSITIVE } from '../../modules/local/seqscreen/taxonomic_assignment/main'
include { SEQSCREEN_FUNCTIONAL_ASSIGNMENT } from '../../modules/local/seqscreen/functional_assignment/main'
include { SEQSCREEN_REPORT_TSV } from '../../modules/local/seqscreen/report_tsv/main'
include { SEQSCREEN_REPORT_HTML } from '../../modules/local/seqscreen/report_html/main'
include { SEQSCREEN_REFERENCE_INFERENCE } from '../../modules/local/seqscreen/reference_inference/main'
include { SEQSCREEN_FORMAT_REPORT } from '../../modules/local/seqscreen/format_report/main'

workflow SEQSCREEN_DSL2 {

    take:
    ch_fasta
    ch_db

    main:
    ch_versions = Channel.empty()
    ch_assets = Channel.value(file("${projectDir}/assets/seqscreen"))

    def mode = (params.seqscreen_mode ?: 'fast').toString().toLowerCase()
    if (!['fast', 'sensitive'].contains(mode)) {
        error("SeqScreen DSL2 currently supports --seqscreen_mode fast or sensitive. Got: ${params.seqscreen_mode}")
    }

    ch_input = ch_fasta.map { record ->
        def values = record instanceof List ? record : [record]
        if (values.size() < 2) {
            error("SEQSCREEN_DSL2 expects [meta, fasta].")
        }

        def meta = values[0]
        def fasta = values[1]
        if (fasta instanceof Collection) {
            if (fasta.size() != 1) {
                error("SeqScreen expects one FASTA per sample '${meta.id}', got ${fasta.size()} files.")
            }
            fasta = fasta[0]
        }
        tuple(meta, fasta)
    }

    SEQSCREEN_VALIDATE_FASTA(ch_input, ch_assets)
    ch_versions = ch_versions.mix(SEQSCREEN_VALIDATE_FASTA.out.versions)

    if (params.seqscreen_window) {
        SEQSCREEN_WINDOW_FASTA(SEQSCREEN_VALIDATE_FASTA.out.fasta, ch_assets)
        ch_seqscreen_fasta = SEQSCREEN_WINDOW_FASTA.out.fasta
        ch_versions = ch_versions.mix(SEQSCREEN_WINDOW_FASTA.out.versions)
    } else {
        ch_seqscreen_fasta = SEQSCREEN_VALIDATE_FASTA.out.fasta
    }

    ch_fasta_key = ch_seqscreen_fasta.map { meta, fasta -> tuple(meta.id, meta, fasta) }

    /*
     * Seq mapper: BSAT/VFDB bowtie2 and rapsearch2 scans.
     */
    SEQSCREEN_BOWTIE2_BSAT(ch_seqscreen_fasta, ch_db, 'bsat')
    SEQSCREEN_BOWTIE2_VFDB(ch_seqscreen_fasta, ch_db, 'vfdb')
    SEQSCREEN_RAPSEARCH2_BSAT(ch_seqscreen_fasta, ch_db, 'bsat')
    SEQSCREEN_RAPSEARCH2_VFDB(ch_seqscreen_fasta, ch_db, 'vfdb')
    ch_versions = ch_versions.mix(
        SEQSCREEN_BOWTIE2_BSAT.out.versions,
        SEQSCREEN_BOWTIE2_VFDB.out.versions,
        SEQSCREEN_RAPSEARCH2_BSAT.out.versions,
        SEQSCREEN_RAPSEARCH2_VFDB.out.versions
    )

    ch_bsat_threats_input = SEQSCREEN_BOWTIE2_BSAT.out.sam
        .map { meta, sam -> tuple(meta.id, meta, sam) }
        .join(SEQSCREEN_RAPSEARCH2_BSAT.out.m8.map { meta, m8 -> tuple(meta.id, meta, m8) }, by: 0)
        .map { sample_id, meta, sam, rap_meta, m8 -> tuple(meta, sam, m8) }

    ch_vfdb_threats_input = SEQSCREEN_BOWTIE2_VFDB.out.sam
        .map { meta, sam -> tuple(meta.id, meta, sam) }
        .join(SEQSCREEN_RAPSEARCH2_VFDB.out.m8.map { meta, m8 -> tuple(meta.id, meta, m8) }, by: 0)
        .map { sample_id, meta, sam, rap_meta, m8 -> tuple(meta, sam, m8) }

    SEQSCREEN_THREATS_BY_BLACKLIST_BSAT(ch_bsat_threats_input, ch_assets, 'blacklist')
    SEQSCREEN_THREATS_BY_BLACKLIST_VFDB(ch_vfdb_threats_input, ch_assets, 'vfdb')
    ch_versions = ch_versions.mix(
        SEQSCREEN_THREATS_BY_BLACKLIST_BSAT.out.versions,
        SEQSCREEN_THREATS_BY_BLACKLIST_VFDB.out.versions
    )

    if (params.seqscreen_hmmscan) {
        SEQSCREEN_HMMSCAN(ch_seqscreen_fasta, ch_db)
        ch_versions = ch_versions.mix(SEQSCREEN_HMMSCAN.out.versions)
    }

    /*
     * Taxonomic and functional annotation.
     */
    if (mode == 'fast') {
        SEQSCREEN_DIAMOND_BLASTX(ch_seqscreen_fasta, ch_db)
        SEQSCREEN_CENTRIFUGE(ch_seqscreen_fasta, ch_db)
        ch_versions = ch_versions.mix(
            SEQSCREEN_DIAMOND_BLASTX.out.versions,
            SEQSCREEN_CENTRIFUGE.out.versions
        )

        ch_tax_input = SEQSCREEN_DIAMOND_BLASTX.out.btab
            .map { meta, btab -> tuple(meta.id, meta, btab) }
            .join(SEQSCREEN_CENTRIFUGE.out.results.map { meta, centrifuge -> tuple(meta.id, meta, centrifuge) }, by: 0)
            .map { sample_id, meta, btab, centrifuge_meta, centrifuge -> tuple(meta, btab, centrifuge) }

        SEQSCREEN_TAXONOMIC_ASSIGNMENT_FAST(ch_tax_input, ch_assets)
        ch_versions = ch_versions.mix(SEQSCREEN_TAXONOMIC_ASSIGNMENT_FAST.out.versions)
        ch_taxonomy = SEQSCREEN_TAXONOMIC_ASSIGNMENT_FAST.out.results
        ch_blastx_btab = SEQSCREEN_DIAMOND_BLASTX.out.functional_btab
        ch_blastx_xml = SEQSCREEN_DIAMOND_BLASTX.out.functional_xml
    } else {
        SEQSCREEN_BLASTN_NT(ch_seqscreen_fasta, ch_db, 'nt')
        SEQSCREEN_BLASTN_MEGARES(ch_seqscreen_fasta, ch_db, 'megares')
        SEQSCREEN_MUMMER(ch_seqscreen_fasta, ch_db)
        SEQSCREEN_BLASTX(ch_seqscreen_fasta, ch_db)
        ch_versions = ch_versions.mix(
            SEQSCREEN_BLASTN_NT.out.versions,
            SEQSCREEN_BLASTN_MEGARES.out.versions,
            SEQSCREEN_MUMMER.out.versions,
            SEQSCREEN_BLASTX.out.versions
        )

        ch_outlier_input = ch_fasta_key
            .join(SEQSCREEN_BLASTN_NT.out.btab.map { meta, btab -> tuple(meta.id, meta, btab) }, by: 0)
            .map { sample_id, meta, fasta, blastn_meta, blastn_btab -> tuple(meta, fasta, blastn_btab) }

        SEQSCREEN_OUTLIER_DETECTION(ch_outlier_input, ch_assets)
        ch_versions = ch_versions.mix(SEQSCREEN_OUTLIER_DETECTION.out.versions)

        ch_tax_input = SEQSCREEN_BLASTX.out.btab
            .map { meta, btab -> tuple(meta.id, meta, btab) }
            .join(SEQSCREEN_OUTLIER_DETECTION.out.clean_btab.map { meta, btab -> tuple(meta.id, meta, btab) }, by: 0)
            .map { sample_id, meta, blastx_btab, outlier_meta, clean_blastn_btab -> tuple(meta, blastx_btab, clean_blastn_btab) }

        SEQSCREEN_TAXONOMIC_ASSIGNMENT_SENSITIVE(ch_tax_input, ch_assets)
        ch_versions = ch_versions.mix(SEQSCREEN_TAXONOMIC_ASSIGNMENT_SENSITIVE.out.versions)
        ch_taxonomy = SEQSCREEN_TAXONOMIC_ASSIGNMENT_SENSITIVE.out.results
        ch_blastx_btab = SEQSCREEN_BLASTX.out.functional_btab
        ch_blastx_xml = SEQSCREEN_BLASTX.out.functional_xml
    }

    ch_functional_input = ch_fasta_key
        .join(ch_blastx_btab.map { meta, btab -> tuple(meta.id, meta, btab) }, by: 0)
        .map { sample_id, meta, fasta, blastx_meta, blastx_btab -> tuple(meta, fasta, blastx_btab) }

    SEQSCREEN_FUNCTIONAL_ASSIGNMENT(ch_functional_input, ch_db, ch_assets)
    ch_versions = ch_versions.mix(SEQSCREEN_FUNCTIONAL_ASSIGNMENT.out.versions)

    /*
     * Reports, reference inference, and final output layout.
     */
    ch_tax_key = ch_taxonomy.map { meta, taxonomy -> tuple(meta.id, meta, taxonomy) }
    ch_func_key = SEQSCREEN_FUNCTIONAL_ASSIGNMENT.out.results.map { meta, functional -> tuple(meta.id, meta, functional) }
    ch_bsat_key = SEQSCREEN_THREATS_BY_BLACKLIST_BSAT.out.threats.map { meta, threats -> tuple(meta.id, meta, threats) }
    ch_vfdb_key = SEQSCREEN_THREATS_BY_BLACKLIST_VFDB.out.threats.map { meta, threats -> tuple(meta.id, meta, threats) }
    ch_xml_key = ch_blastx_xml.map { meta, xml -> tuple(meta.id, meta, xml) }

    ch_report_input = ch_tax_key
        .join(ch_func_key, by: 0)
        .join(ch_bsat_key, by: 0)
        .join(ch_vfdb_key, by: 0)
        .join(ch_fasta_key, by: 0)
        .map { sample_id, meta, taxonomy, func_meta, functional, bsat_meta, bsat, vfdb_meta, vfdb, fasta_meta, fasta ->
            tuple(meta, fasta, taxonomy, functional, bsat, vfdb)
        }

    SEQSCREEN_REPORT_TSV(ch_report_input, ch_db, ch_assets, mode)
    ch_versions = ch_versions.mix(SEQSCREEN_REPORT_TSV.out.versions)

    ch_report_key = SEQSCREEN_REPORT_TSV.out.report.map { meta, report -> tuple(meta.id, meta, report) }

    ch_html_input = ch_report_key
        .join(ch_xml_key, by: 0)
        .join(ch_fasta_key, by: 0)
        .map { sample_id, meta, report, xml_meta, blastx_xml, fasta_meta, fasta ->
            tuple(meta, fasta, report, blastx_xml)
        }

    def rflag = mode
    SEQSCREEN_REPORT_HTML(ch_html_input, ch_db, ch_assets, mode, rflag)
    ch_versions = ch_versions.mix(SEQSCREEN_REPORT_HTML.out.versions)

    ch_reference_input = ch_report_key
        .join(ch_fasta_key, by: 0)
        .map { sample_id, meta, report, fasta_meta, fasta -> tuple(meta, fasta, report) }

    SEQSCREEN_REFERENCE_INFERENCE(ch_reference_input, ch_db, ch_assets)
    ch_versions = ch_versions.mix(SEQSCREEN_REFERENCE_INFERENCE.out.versions)

    ch_html_key = SEQSCREEN_REPORT_HTML.out.html.map { meta, html -> tuple(meta.id, meta, html) }
    ch_reference_key = SEQSCREEN_REFERENCE_INFERENCE.out.outdir.map { meta, outdir -> tuple(meta.id, meta, outdir) }

    ch_format_input = ch_report_key
        .join(ch_html_key, by: 0)
        .join(ch_tax_key, by: 0)
        .join(ch_func_key, by: 0)
        .join(ch_reference_key, by: 0)
        .map { sample_id, meta, report, html_meta, html, tax_meta, taxonomy, func_meta, functional, ref_meta, reference_inference ->
            tuple(meta, report, html, taxonomy, functional, reference_inference)
        }

    SEQSCREEN_FORMAT_REPORT(ch_format_input, ch_db, ch_assets, mode)
    ch_versions = ch_versions.mix(SEQSCREEN_FORMAT_REPORT.out.versions)

    emit:
    output_dir         = SEQSCREEN_FORMAT_REPORT.out.outdir
    report             = SEQSCREEN_FORMAT_REPORT.out.report
    html               = SEQSCREEN_FORMAT_REPORT.out.html
    taxonomic_results  = SEQSCREEN_FORMAT_REPORT.out.taxonomic_results
    functional_results = SEQSCREEN_FORMAT_REPORT.out.functional_results
    versions           = ch_versions
}
