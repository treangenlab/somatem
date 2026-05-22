#!/usr/bin/env nextflow

// Metagenomics analysis subworkflow: taxonomic profiling, de novo assembly, mapping, binning, quality assessment, and annotation

  
// Include nf-core modules
include { FLYE }                    from '../../modules/nf-core/flye/main'
include { MINIMAP2_INDEX }          from '../../modules/nf-core/minimap2/index/main'
include { MINIMAP2_ALIGN }          from '../../modules/nf-core/minimap2/align/main'
include { SAMTOOLS_SORT }           from '../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_COVERAGE }       from '../../modules/nf-core/samtools/coverage/main'
include { SEMIBIN_SINGLEEASYBIN }   from '../../modules/nf-core/semibin/singleeasybin/main'
include { CHECKM2_PREDICT }         from '../../modules/nf-core/checkm2/predict/main'
include { CHECKM2_PARSE }           from '../../modules/local/checkm2/parse/main'
include { BAKTA_BAKTA }             from '../../modules/local/bakta/bakta/main'
include { SINGLEM_PIPE }            from '../../modules/local/singlem/pipe/main'
include { SINGLEM_PIPE as SINGLEM_PIPE_BINS } from '../../modules/local/singlem/pipe/main'
include { SINGLEM_APPRAISE }        from '../../modules/local/singlem/appraise/main'
include { TAXBURST }                from '../../modules/local/taxburst/main'
include { PIGEON }                  from '../../modules/local/pigeon/main'
include { PIGEON_ITERATIVE_BINNING } from '../../modules/local/pigeon/iterate/main'
include { AGB }                     from '../../modules/local/agb/main'
include { SOMATEM_SUMMARY_REPORT as ASSEMBLY_MAGS_SUMMARY_REPORT } from '../../modules/local/somatem_summary_report/main.nf'


workflow ASSEMBLY_MAGS {

    take:
    reads               // channel: [ val(meta), path(reads) ]
    
    ch_checkm2_db       // channel: path(db)
    ch_bakta_db         // channel: path(db)
    ch_singlem_db // channel: path(metapackage)

    main:
    ch_versions = channel.empty()

    // Check input channel (debug)
    // reads.view { meta, file -> "Input reads: ${meta.id} -> ${file}" } // debug

    // Taxonomic profiling with SingleM on raw reads
    SINGLEM_PIPE(reads, ch_singlem_db, 'reads')
    ch_versions = ch_versions.mix(SINGLEM_PIPE.out.versions)

    // Interactive taxonomic visualization with TaxBurst
    TAXBURST(SINGLEM_PIPE.out.taxonomic_profile, 'SingleM')
    ch_versions = ch_versions.mix(TAXBURST.out.versions)
    TAXBURST.out.html.view { meta, _html -> "✓ Taxonomic profiling and interactive visualization completed for ${meta.id}" } // log

    // Assembly with Flye
    FLYE(reads, params.flye_mode)
    ch_versions = ch_versions.mix(FLYE.out.versions)
    FLYE.out.fasta.view { meta, _fasta -> "✓ Assembly completed for ${meta.id}" } // log

    // visualize assembly graph with agb
    AGB(FLYE.out.gfa, FLYE.out.gv, FLYE.out.txt)
    AGB.out.assembly_graph.view { meta, _graph -> "✓ Assembly graph visualization completed for ${meta.id}" } // log

    // Create minimap2 index
    MINIMAP2_INDEX(FLYE.out.fasta)
    ch_versions = ch_versions.mix(MINIMAP2_INDEX.out.versions)
    MINIMAP2_INDEX.out.index.view { meta, _index -> "✓ Minimap2 index created for ${meta.id}" } // log  

    // Map reads back to assembly
    MINIMAP2_ALIGN(
        reads,
        MINIMAP2_INDEX.out.index,
        true,
        'bai',
        false,
        false
    )
    ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions)

    // Sort BAM files
    SAMTOOLS_SORT(MINIMAP2_ALIGN.out.bam, FLYE.out.fasta, 'bai')
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions)

    // Calculate coverage
    ch_bam_for_coverage = SAMTOOLS_SORT.out.bam.map { meta, bam ->
        [meta, bam, []]
    }
    
    ch_fasta_for_coverage = FLYE.out.fasta
    
    ch_fqi_for_coverage = SAMTOOLS_SORT.out.bam.map { meta, _bam ->
        [meta, file('OPTIONAL_FILE')]
    }
    
    SAMTOOLS_COVERAGE(
        ch_bam_for_coverage,
        ch_fasta_for_coverage,
        ch_fqi_for_coverage
    )
    ch_versions = ch_versions.mix(SAMTOOLS_COVERAGE.out.versions)
    SAMTOOLS_COVERAGE.out.coverage.view { meta, _coverage -> "✓ Coverage calculated for ${meta.id}" } // log

    // Binning with Pigeon-controlled iterative binning by default, or plain SemiBin2 when disabled.
    use_iterative_binning = params.mag_iterative_binning_enabled instanceof Boolean
        ? params.mag_iterative_binning_enabled
        : params.mag_iterative_binning_enabled.toString().toBoolean()

    ch_asm_bam = FLYE.out.fasta.join(SAMTOOLS_SORT.out.bam, by: [0]) // join by sample ID (meta.id)

    ch_bins_for_downstream = channel.empty()
    ch_bins_csv = channel.empty()
    ch_bins_tsv = channel.empty()
    ch_bins_model = channel.empty()
    ch_iterative_dastool_bins = channel.empty()
    ch_iterative_selection_summary = channel.empty()
    ch_iterative_manifest = channel.empty()
    ch_iterative_selected = channel.empty()
    ch_iterative_trajectory = channel.empty()
    ch_iterative_summary = channel.empty()
    ch_iterative_report = channel.empty()
    ch_iterative_command_log = channel.empty()

    if (!use_iterative_binning) {
        SEMIBIN_SINGLEEASYBIN(ch_asm_bam)
        ch_versions = ch_versions.mix(SEMIBIN_SINGLEEASYBIN.out.versions)
        SEMIBIN_SINGLEEASYBIN.out.csv.view { meta, _csv -> "✓ Binning completed for ${meta.id}" } // log
        ch_bins_for_downstream = SEMIBIN_SINGLEEASYBIN.out.output_fasta
        ch_bins_csv = SEMIBIN_SINGLEEASYBIN.out.csv
        ch_bins_tsv = SEMIBIN_SINGLEEASYBIN.out.tsv
        ch_bins_model = SEMIBIN_SINGLEEASYBIN.out.model
    } else {
        ch_iterative_input = FLYE.out.gfa
            .join(ch_asm_bam, by: [0])
            .map { meta, gfa, assembly, bam ->
                [meta, gfa, assembly, bam]
            }

        PIGEON_ITERATIVE_BINNING(ch_iterative_input)
        ch_versions = ch_versions.mix(PIGEON_ITERATIVE_BINNING.out.versions)

        // With DAS Tool consensus enabled, send the consensus bins directly to CheckM2.
        // Otherwise, use final_bins from the selected binner iterations.
        ch_bins_for_downstream = params.mag_iterative_consensus_tool == "dastool"
            ? PIGEON_ITERATIVE_BINNING.out.dastool_bins
            : PIGEON_ITERATIVE_BINNING.out.final_bins
        ch_iterative_dastool_bins = PIGEON_ITERATIVE_BINNING.out.dastool_bins
        ch_iterative_selection_summary = PIGEON_ITERATIVE_BINNING.out.selection_summary
        ch_iterative_manifest = PIGEON_ITERATIVE_BINNING.out.candidate_manifest
        ch_iterative_selected = PIGEON_ITERATIVE_BINNING.out.selected_tsv
        ch_iterative_trajectory = PIGEON_ITERATIVE_BINNING.out.trajectory_tsv
        ch_iterative_summary = PIGEON_ITERATIVE_BINNING.out.iterative_summary
        ch_iterative_report = PIGEON_ITERATIVE_BINNING.out.iterative_report
        ch_iterative_command_log = PIGEON_ITERATIVE_BINNING.out.command_log
        PIGEON_ITERATIVE_BINNING.out.iterative_summary.view { meta, _summary -> "✓ Optional iterative binning completed for ${meta.id}" } // log
    }

    // Quality assessment with CheckM2
    CHECKM2_PREDICT(ch_bins_for_downstream, ch_checkm2_db)
    ch_versions = ch_versions.mix(CHECKM2_PREDICT.out.versions)
    CHECKM2_PREDICT.out.checkm2_tsv.view { meta, _tsv -> "✓ Quality assessment completed for ${meta.id}" } // log
    
    // Parse CheckM2 results to get completeness information
    CHECKM2_PARSE(CHECKM2_PREDICT.out.checkm2_tsv, ch_bins_for_downstream)
    ch_versions = ch_versions.mix(CHECKM2_PARSE.out.versions)

    // Run SingleM pipe on bins - FIXED: Add suffix to avoid filename collision
    SINGLEM_PIPE_BINS(ch_bins_for_downstream, ch_singlem_db, 'genome')
    ch_versions = ch_versions.mix(SINGLEM_PIPE_BINS.out.versions)

    
    // PIGEON ANALYSIS: compare k-mer composition from unitigs, contigs and bins
    ch_pigeon_input = FLYE.out.gfa // join gfa, assembly, and bins to prepare pigeon input
        .join(FLYE.out.fasta, by: [0])
        .join(ch_bins_for_downstream, by: [0])

    PIGEON(ch_pigeon_input)
    ch_versions = ch_versions.mix(PIGEON.out.versions)
    PIGEON.out.report.view { meta, _report -> "✓ Pigeon post-hoc analysis completed for ${meta.id}" } // log


    // SingleM appraise - simplified input preparation
    ch_metagenome_otu = SINGLEM_PIPE.out.otu_table
        .map { meta, otu -> 
            def clean_meta = [id: meta.id]
            return [clean_meta, otu]
        }
    
    ch_bins_otu = SINGLEM_PIPE_BINS.out.otu_table
        .map { meta, otu -> 
            // Extract original sample ID by removing the "_bins" suffix
            def original_id = meta.id.replaceAll(/_bins$/, '')
            def clean_meta = [id: original_id]
            return [clean_meta, otu]
        }

    ch_appraise_input = ch_metagenome_otu
        .cross(ch_bins_otu) { it[0].id }
        .map { metagenome_tuple, bins_tuple ->
            def meta = metagenome_tuple[0]
            def metOtu = metagenome_tuple[1]
            def binOtu = bins_tuple[1]
            
            def new_meta = [id: "${meta.id}_appraisal"]
            
            // Handle potential ArrayList/Collection objects by extracting the actual file
            def actualMetOtu = (metOtu instanceof Collection && metOtu.size() > 0) ? metOtu[0] : metOtu
            def actualBinOtu = (binOtu instanceof Collection && binOtu.size() > 0) ? binOtu[0] : binOtu
            
            [new_meta, actualMetOtu, actualBinOtu, []]
        }

    SINGLEM_APPRAISE(ch_appraise_input, ch_singlem_db)
    ch_versions = ch_versions.mix(SINGLEM_APPRAISE.out.versions)
    SINGLEM_APPRAISE.out.summary.view { meta, _summary -> "✓ SingleM appraise analysis completed for ${meta.id}" } // log

    // Completeness-based Bakta annotation
    ch_bins_for_annotation = ch_bins_for_downstream
        .transpose()
        .map { meta, bin ->
            def new_meta = meta.clone()
            
            // Get the actual bin filename and clean it
            def bin_filename = bin.name
            def clean_bin_name = bin_filename.replaceAll(/\.(fa|fasta|fna)(\.gz)?$/, '')
            
            // Set the meta.id to the clean bin name and preserve sample_id for folder organization
            new_meta.id = clean_bin_name
            new_meta.sample_id = meta.id  // Keep original sample ID for folder structure
            
            // Extract the SemiBin part for completeness lookup
            def semibin_match = bin_filename =~ /.*?(SemiBin_\d+)\.fa/
            def semibin_name = semibin_match ? semibin_match[0][1] : clean_bin_name
            new_meta.bin_name = semibin_name
            
            [new_meta, bin]
        }

    // Add completeness information to bin metadata
    ch_completeness_with_meta = CHECKM2_PARSE.out.completeness_map
        .map { meta, csv_file ->
            // Parse the CSV file to create a map but keep the metadata
            def completeness_map = [:]
            csv_file.text.split('\n').drop(1).each { line ->
                if (line.trim()) {
                    def parts = line.split(',')
                    if (parts.size() >= 2) {
                        completeness_map[parts[0]] = Float.parseFloat(parts[1])
                    }
                }
            }
            return [meta, completeness_map]
        }

    // Use sample_id for joining completeness data with bins
    ch_bins_with_completeness = ch_bins_for_annotation
        .map { meta, bin -> 
            def join_key = meta.sample_id
            [join_key, meta, bin]
        }
        .combine(
            ch_completeness_with_meta.map { meta, completeness_map -> 
                [meta.id, completeness_map] 
            }, 
            by: 0
        )
        .map { _join_key, meta, bin, completeness_map ->
            def new_meta = meta.clone()
            def bin_name = meta.bin_name
            
            // Look up completeness for this bin
            def completeness = completeness_map[bin_name]
            if (completeness != null) {
                new_meta.completeness = completeness
            } else {
                new_meta.completeness = null
            }
            
            [new_meta, bin]
        }
        .filter { meta, _bin ->
            // Only process bins that meet the completeness threshold
            return meta.completeness != null && meta.completeness >= params.checkm2_completeness_threshold
        }

    // Only run Bakta on high-quality bins (≥checkm2_completeness_threshold)
    BAKTA_BAKTA(ch_bins_with_completeness, ch_bakta_db, [], [])
    ch_versions = ch_versions.mix(BAKTA_BAKTA.out.versions)

    // Show which bins got annotation (log)
    BAKTA_BAKTA.out.embl.view { meta, _embl -> 
        def completeness = meta.completeness ?: "unknown"
        "✓ Bakta annotation completed for ${meta.id} (${completeness}% complete)"
    }

    // Show which bins got annotation (log)
    BAKTA_BAKTA.out.embl
        .map { meta, _embl -> meta.completeness }
        .filter { it != null && it >= params.checkm2_completeness_threshold }
        .count()
        .view { count -> "✓ Generated annotations for ${count} high-quality bins (≥${params.checkm2_completeness_threshold}% complete)" }

    ch_versions_for_report = ch_versions
    ch_assembly_report_files = FLYE.out.txt
        .mix(
            FLYE.out.log,
            AGB.out.assembly_graph,
            SAMTOOLS_COVERAGE.out.coverage,
            ch_bins_csv,
            ch_bins_tsv,
            ch_iterative_manifest,
            ch_iterative_selected,
            ch_iterative_trajectory,
            ch_iterative_summary,
            ch_iterative_report,
            CHECKM2_PREDICT.out.checkm2_tsv,
            CHECKM2_PARSE.out.completeness_map,
            SINGLEM_PIPE.out.taxonomic_profile,
            SINGLEM_PIPE_BINS.out.otu_table,
            TAXBURST.out.html,
            PIGEON.out.report,
            PIGEON.out.metrics,
            SINGLEM_APPRAISE.out.summary,
            SINGLEM_APPRAISE.out.plot,
            BAKTA_BAKTA.out.tsv,
            BAKTA_BAKTA.out.json,
            BAKTA_BAKTA.out.png,
            ch_versions_for_report
        )

    ASSEMBLY_MAGS_SUMMARY_REPORT(
        'assembly_mags',
        'Assembly and MAG recovery',
        Channel.fromPath(params.input),
        ch_assembly_report_files.flatMap { item ->
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
    ch_versions = ch_versions.mix(ASSEMBLY_MAGS_SUMMARY_REPORT.out.versions)


    emit:
    // Original outputs
    assembly        = FLYE.out.fasta
    assembly_gfa    = FLYE.out.gfa
    assembly_log    = FLYE.out.log
    assembly_graph  = AGB.out.assembly_graph
    minimap_index   = MINIMAP2_INDEX.out.index
    bam_aligned     = MINIMAP2_ALIGN.out.bam
    bam_aligned_index = MINIMAP2_ALIGN.out.index
    bam_sorted      = SAMTOOLS_SORT.out.bam
    bam_sorted_index = SAMTOOLS_SORT.out.bai
    coverage        = SAMTOOLS_COVERAGE.out.coverage
    bins            = ch_bins_for_downstream
    bins_csv        = ch_bins_csv
    bins_tsv        = ch_bins_tsv
    bins_model      = ch_bins_model
    iterative_dastool_bins = ch_iterative_dastool_bins
    iterative_selection_summary = ch_iterative_selection_summary
    iterative_manifest = ch_iterative_manifest
    iterative_selected = ch_iterative_selected
    iterative_trajectory = ch_iterative_trajectory
    iterative_summary = ch_iterative_summary
    iterative_report = ch_iterative_report
    iterative_command_log = ch_iterative_command_log
    checkm2_report  = CHECKM2_PREDICT.out.checkm2_tsv
    checkm2_output  = CHECKM2_PREDICT.out.checkm2_output
    completeness_map = CHECKM2_PARSE.out.completeness_map
    
    // Pigeon outputs
    pigeon_html     = PIGEON.out.report
    pigeon_metrics  = PIGEON.out.metrics

    // Bakta annotation outputs
    bakta_embl              = BAKTA_BAKTA.out.embl
    bakta_faa               = BAKTA_BAKTA.out.faa
    bakta_ffn               = BAKTA_BAKTA.out.ffn
    bakta_fna               = BAKTA_BAKTA.out.fna
    bakta_gbff              = BAKTA_BAKTA.out.gbff
    bakta_gff               = BAKTA_BAKTA.out.gff
    bakta_hypotheticals_tsv = BAKTA_BAKTA.out.hypotheticals_tsv
    bakta_hypotheticals_faa = BAKTA_BAKTA.out.hypotheticals_faa
    bakta_tsv               = BAKTA_BAKTA.out.tsv
    bakta_txt               = BAKTA_BAKTA.out.txt
    bakta_inference_tsv     = BAKTA_BAKTA.out.inference_tsv
    bakta_png               = BAKTA_BAKTA.out.png
    bakta_svg               = BAKTA_BAKTA.out.svg
    bakta_json              = BAKTA_BAKTA.out.json
    
    // Taxonomic profiling outputs
    singlem_profile = SINGLEM_PIPE.out.taxonomic_profile
    singlem_otu     = SINGLEM_PIPE.out.otu_table
    singlem_bins_otu = SINGLEM_PIPE_BINS.out.otu_table
    taxburst_html   = TAXBURST.out.html
    
    // SingleM appraise outputs
    appraise_summary = SINGLEM_APPRAISE.out.summary
    appraise_plot = SINGLEM_APPRAISE.out.plot
    appraise_binned_otu = SINGLEM_APPRAISE.out.binned_otu_table
    appraise_unbinned_otu = SINGLEM_APPRAISE.out.unbinned_otu_table
    appraise_assembled_otu = SINGLEM_APPRAISE.out.assembled_otu_table
    appraise_unaccounted_otu = SINGLEM_APPRAISE.out.unaccounted_otu_table
    
    summary_report  = ASSEMBLY_MAGS_SUMMARY_REPORT.out.html
    versions        = ch_versions
}
