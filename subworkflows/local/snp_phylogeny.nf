#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Reference-based whole-genome SNP phylogeny
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Nextflow-native route for assembled isolate genomes:
      genomes + reference -> minibwa PAF evidence + Parsnp core SNP alignment -> RAxML-NG
----------------------------------------------------------------------------------------
*/

include { MINIBWA_INDEX } from '../../modules/local/minibwa/index/main.nf'
include { MINIBWA_MAP } from '../../modules/local/minibwa/map/main.nf'
include { PARSNP } from '../../modules/local/parsnp/parsnp/main.nf'
include { RAXMLNG_TREE } from '../../modules/local/raxmlng/tree/main.nf'
include { SOMATEM_SUMMARY_REPORT as SNP_PHYLOGENY_SUMMARY_REPORT } from '../../modules/local/somatem_summary_report/main.nf'

workflow SNP_PHYLOGENY {

    take:
    ch_genomes

    main:
    ch_versions = channel.empty()

    if (!params.reference) {
        error("SNP phylogeny requires --reference with a reference genome FASTA.")
    }

    ch_reference = Channel.value(file(params.reference, checkIfExists: true))

    MINIBWA_INDEX(ch_reference)
    ch_versions = ch_versions.mix(MINIBWA_INDEX.out.versions)

    MINIBWA_MAP(ch_genomes, MINIBWA_INDEX.out.index)
    ch_versions = ch_versions.mix(MINIBWA_MAP.out.versions)

    ch_query_genomes = ch_genomes.map { meta, genome -> genome }.collect()

    PARSNP(ch_reference, ch_query_genomes)
    ch_versions = ch_versions.mix(PARSNP.out.versions)

    RAXMLNG_TREE(PARSNP.out.core_snps_fasta)
    ch_versions = ch_versions.mix(RAXMLNG_TREE.out.versions)

    ch_snp_phylo_outputs = MINIBWA_MAP.out.paf
        .mix(
            PARSNP.out.outdir,
            PARSNP.out.gingr_archive,
            PARSNP.out.xmfa,
            PARSNP.out.tree,
            PARSNP.out.core_snps_fasta,
            RAXMLNG_TREE.out.tree,
            RAXMLNG_TREE.out.raxmlng_dir
        )

    ch_versions_for_report = ch_versions
    SNP_PHYLOGENY_SUMMARY_REPORT(
        'snp_phylogeny',
        'Reference-based SNP phylogeny',
        Channel.fromPath(params.input),
        ch_snp_phylo_outputs.mix(ch_versions_for_report).flatMap { item ->
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
    ch_versions = ch_versions.mix(SNP_PHYLOGENY_SUMMARY_REPORT.out.versions)

    emit:
    alignments          = MINIBWA_MAP.out.paf
    parsnp_outdir       = PARSNP.out.outdir
    parsnp_ggr          = PARSNP.out.gingr_archive
    parsnp_xmfa         = PARSNP.out.xmfa
    parsnp_tree         = PARSNP.out.tree
    core_snps_fasta     = PARSNP.out.core_snps_fasta
    tree                = RAXMLNG_TREE.out.tree
    raxmlng_dir         = RAXMLNG_TREE.out.raxmlng_dir
    summary_report      = SNP_PHYLOGENY_SUMMARY_REPORT.out.html
    versions            = ch_versions
}
