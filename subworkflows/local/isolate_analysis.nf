#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Long-read-first bacterial isolate analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Initial isolate workflow route:

      ONT/PacBio long reads
        -> Autocycler long-read subsampling
        -> multiple Flye assemblies
        -> Autocycler consensus assembly
        -> final long-read-first assembly

      Optional hybrid polishing when params.hybrid_assembly is true:
        -> BWA short-read alignments for Polypolish
        -> Polypolish
        -> Pypolca --careful
        -> final polished hybrid assembly

    Expected input channel:
      tuple val(meta), path(long_reads)
      OR
      tuple val(meta), path(long_reads), path(short_reads_1), path(short_reads_2),
            val(expected_genome_size), val(species)
      For long-read-only assembly, short_reads_1 and short_reads_2 may be empty
      placeholders because they are only consumed when params.hybrid_assembly is true.

    Expected meta:
      meta.id is required.
      Optional meta.shortread_depth / meta.short_read_depth / meta.illumina_depth can be
      used to skip Pypolca below params.min_shortread_depth_for_pypolca.
----------------------------------------------------------------------------------------
*/

include { AUTOCYCLER_SUBSAMPLE }   from '../../modules/local/autocycler/subsample/main'
include { FLYE_AUTOCYCLER }        from '../../modules/local/flye/autocycler/main'
include { AUTOCYCLER_CONSENSUS }   from '../../modules/local/autocycler/consensus/main'
include { BWA_MEM_FOR_POLYPOLISH } from '../../modules/local/bwa/mem_for_polypolish/main'
include { POLYPOLISH_FILTER }      from '../../modules/local/polypolish/filter/main'
include { POLYPOLISH_POLISH }      from '../../modules/local/polypolish/polish/main'
include { PYPOLCA }                from '../../modules/local/pypolca/main'
include { FASTA_FINALIZE }         from '../../modules/local/fasta/finalize/main'
include { KRAKEN2_KRAKEN2 }        from '../../modules/nf-core/kraken2/kraken2/main'
include { BAKTA_BAKTA }            from '../../modules/local/bakta/bakta/main'
include { BTYPER3 }                from '../../modules/local/btyper3/main'
include { SOMATEM_SUMMARY_REPORT as ISOLATE_ANALYSIS_SUMMARY_REPORT } from '../../modules/local/somatem_summary_report/main.nf'

workflow ISOLATE_ANALYSIS {

    take:
    ch_reads
    ch_bakta_db
    ch_kraken2_db
    // tuple val(meta), path(long_reads), path(short_reads_1), path(short_reads_2),
    //       val(expected_genome_size), val(species)

    main:
    ch_versions = channel.empty()

    def hybrid_assembly = params.hybrid_assembly instanceof Boolean
        ? params.hybrid_assembly
        : (params.hybrid_assembly == null ? false : params.hybrid_assembly.toString().toBoolean())

    def run_btyper3 = params.run_btyper3 instanceof Boolean
        ? params.run_btyper3
        : (params.run_btyper3 == null ? true : params.run_btyper3.toString().toBoolean())

    ch_isolate_input = ch_reads.map { record ->
        def values = record instanceof List ? record : [record]
        if (values.size() == 2) {
            return tuple(values[0], values[1], '', '', '', '')
        }
        if (values.size() >= 6) {
            return tuple(values[0], values[1], values[2] ?: '', values[3] ?: '', values[4] ?: '', values[5] ?: '')
        }
        error("ISOLATE_ANALYSIS expects [meta, long_reads] or [meta, long_reads, short_reads_1, short_reads_2, genome_size, species].")
    }

    ch_kraken_reads = ch_isolate_input.map { meta, long_reads, short_reads_1, short_reads_2, genome_size, species ->
        def kraken_meta = meta.single_end == null ? meta + [single_end: true] : meta
        tuple(kraken_meta, long_reads)
    }

    KRAKEN2_KRAKEN2(
        ch_kraken_reads,
        ch_kraken2_db,
        params.kraken2_save_output_fastqs,
        params.kraken2_save_reads_assignment
    )
    ch_versions = ch_versions.mix(KRAKEN2_KRAKEN2.out.versions)

    ch_autocycler_input = ch_isolate_input.map { meta, long_reads, short_reads_1, short_reads_2, genome_size, species ->
        tuple(meta, long_reads, genome_size ?: params.autocycler_genome_size ?: '', species ?: '')
    }

    AUTOCYCLER_SUBSAMPLE(ch_autocycler_input)
    ch_versions = ch_versions.mix(AUTOCYCLER_SUBSAMPLE.out.versions)

    ch_subsampled_reads = AUTOCYCLER_SUBSAMPLE.out.reads.flatMap { meta, subsample_reads, genome_size_file, species ->
        def reads = subsample_reads instanceof List ? subsample_reads : [subsample_reads]
        reads.withIndex().collect { read, idx ->
            def subsample_meta = meta + [subsample_id: String.format('subsample_%02d', idx + 1)]
            tuple(subsample_meta, read, genome_size_file, species)
        }
    }

    FLYE_AUTOCYCLER(ch_subsampled_reads)
    ch_versions = ch_versions.mix(FLYE_AUTOCYCLER.out.versions)

    ch_grouped_assemblies = FLYE_AUTOCYCLER.out.assembly
        .map { meta, assembly -> tuple(meta.id, meta, assembly) }
        .groupTuple(by: 0)
        .map { sample_id, metas, assemblies ->
            def base_meta = metas[0].findAll { key, value -> !['subsample_id', 'seed'].contains(key) }
            tuple(base_meta, assemblies)
        }

    AUTOCYCLER_CONSENSUS(ch_grouped_assemblies)
    ch_versions = ch_versions.mix(AUTOCYCLER_CONSENSUS.out.versions)

    def run_polypolish = hybrid_assembly && (params.run_polypolish instanceof Boolean
        ? params.run_polypolish
        : (params.run_polypolish == null ? true : params.run_polypolish.toString().toBoolean()))

    def run_pypolca = hybrid_assembly && (params.run_pypolca instanceof Boolean
        ? params.run_pypolca
        : (params.run_pypolca == null ? true : params.run_pypolca.toString().toBoolean()))

    if (hybrid_assembly) {
        ch_short_reads = ch_isolate_input.map { meta, long_reads, short_reads_1, short_reads_2, genome_size, species ->
            if (!short_reads_1 || !short_reads_2) {
                error("Hybrid isolate assembly requires paired short reads for sample '${meta.id}'.")
            }
            tuple(meta.id, meta, short_reads_1, short_reads_2, genome_size ?: '', species ?: '')
        }

        ch_autocycler_with_short_reads = AUTOCYCLER_CONSENSUS.out.assembly
            .map { meta, assembly -> tuple(meta.id, meta, assembly) }
            .join(ch_short_reads, by: 0)
            .map { sample_id, meta, assembly, read_meta, short_reads_1, short_reads_2, genome_size, species ->
                tuple(meta, assembly, short_reads_1, short_reads_2)
            }
    } else {
        ch_autocycler_with_short_reads = channel.empty()
        ch_short_reads = channel.empty()
    }

    if (hybrid_assembly && !run_polypolish && run_pypolca) {
        ch_polished_with_short_reads = AUTOCYCLER_CONSENSUS.out.assembly
            .map { meta, assembly -> tuple(meta.id, meta, assembly) }
            .join(ch_short_reads, by: 0)
            .map { sample_id, meta, assembly, read_meta, short_reads_1, short_reads_2, genome_size, species ->
                tuple(meta, assembly, short_reads_1, short_reads_2)
            }
    } else {
        ch_polished_with_short_reads = channel.empty()
    }

    if (!hybrid_assembly) {
        ch_polished_assembly = AUTOCYCLER_CONSENSUS.out.assembly
    } else if (!run_polypolish) {
        ch_polished_assembly = AUTOCYCLER_CONSENSUS.out.assembly
    }

    if (hybrid_assembly && run_polypolish) {
        BWA_MEM_FOR_POLYPOLISH(ch_autocycler_with_short_reads)
        ch_versions = ch_versions.mix(BWA_MEM_FOR_POLYPOLISH.out.versions)

        POLYPOLISH_FILTER(BWA_MEM_FOR_POLYPOLISH.out.alignments)
        ch_versions = ch_versions.mix(POLYPOLISH_FILTER.out.versions)

        POLYPOLISH_POLISH(POLYPOLISH_FILTER.out.alignments)
        ch_versions = ch_versions.mix(POLYPOLISH_POLISH.out.versions)
        ch_polished_assembly = POLYPOLISH_POLISH.out.assembly
    }

    if (hybrid_assembly && run_pypolca) {
        def min_shortread_depth = (params.min_shortread_depth_for_pypolca ?: 25) as double

        if (run_polypolish) {
            ch_polished_with_short_reads = ch_polished_assembly
                .map { meta, assembly -> tuple(meta.id, meta, assembly) }
                .join(ch_short_reads, by: 0)
                .map { sample_id, meta, assembly, read_meta, short_reads_1, short_reads_2, genome_size, species ->
                    tuple(meta, assembly, short_reads_1, short_reads_2)
                }
        }

        ch_pypolca_input = ch_polished_with_short_reads.filter { meta, assembly, short_reads_1, short_reads_2 ->
            def depth = meta.shortread_depth ?: meta.short_read_depth ?: meta.illumina_depth
            depth == null || (depth as double) >= min_shortread_depth
        }

        ch_skip_pypolca = ch_polished_with_short_reads
            .filter { meta, assembly, short_reads_1, short_reads_2 ->
                def depth = meta.shortread_depth ?: meta.short_read_depth ?: meta.illumina_depth
                depth != null && (depth as double) < min_shortread_depth
            }
            .map { meta, assembly, short_reads_1, short_reads_2 -> tuple(meta, assembly) }

        PYPOLCA(ch_pypolca_input)
        ch_versions = ch_versions.mix(PYPOLCA.out.versions)
        ch_final_candidate_assembly = PYPOLCA.out.assembly.mix(ch_skip_pypolca)
    } else {
        ch_final_candidate_assembly = ch_polished_assembly
    }

    FASTA_FINALIZE(ch_final_candidate_assembly)
    ch_versions = ch_versions.mix(FASTA_FINALIZE.out.versions)

    ch_bakta_input = FASTA_FINALIZE.out.assembly.map { meta, assembly ->
        def bakta_meta = meta.sample_id == null ? meta + [sample_id: meta.id] : meta
        tuple(bakta_meta, assembly)
    }

    BAKTA_BAKTA(ch_bakta_input, ch_bakta_db, [], [])
    ch_versions = ch_versions.mix(BAKTA_BAKTA.out.versions)

    if (run_btyper3) {
        BTYPER3(BAKTA_BAKTA.out.fna)
        ch_versions = ch_versions.mix(BTYPER3.out.versions)
        ch_btyper3_outdir = BTYPER3.out.outdir
    } else {
        ch_btyper3_outdir = channel.empty()
    }

    ch_versions_for_report = ch_versions
    ch_isolate_report_files = KRAKEN2_KRAKEN2.out.report
        .mix(
            FASTA_FINALIZE.out.assembly,
            AUTOCYCLER_CONSENSUS.out.assembly,
            AUTOCYCLER_CONSENSUS.out.graph,
            BAKTA_BAKTA.out.tsv,
            BAKTA_BAKTA.out.json,
            BAKTA_BAKTA.out.png,
            ch_btyper3_outdir,
            ch_versions_for_report
        )

    ISOLATE_ANALYSIS_SUMMARY_REPORT(
        'isolate_analysis',
        'Isolate assembly, classification, annotation, and typing',
        Channel.fromPath(params.input),
        ch_isolate_report_files.flatMap { item ->
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
    ch_versions = ch_versions.mix(ISOLATE_ANALYSIS_SUMMARY_REPORT.out.versions)

    emit:
    assembly             = FASTA_FINALIZE.out.assembly
    autocycler_assembly  = AUTOCYCLER_CONSENSUS.out.assembly
    autocycler_graph     = AUTOCYCLER_CONSENSUS.out.graph
    autocycler_outdir    = AUTOCYCLER_CONSENSUS.out.outdir
    kraken2_report       = KRAKEN2_KRAKEN2.out.report
    bakta_tsv            = BAKTA_BAKTA.out.tsv
    bakta_json           = BAKTA_BAKTA.out.json
    bakta_gff            = BAKTA_BAKTA.out.gff
    btyper3_outdir       = ch_btyper3_outdir
    summary_report       = ISOLATE_ANALYSIS_SUMMARY_REPORT.out.html
    versions             = ch_versions
}
