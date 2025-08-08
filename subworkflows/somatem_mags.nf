#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Metagenomics analysis subworkflow: taxonomic profiling first, then assembly, mapping, binning, quality assessment, and annotation

// Include nf-core modules
include { FLYE }                    from '../modules/nf-core/flye/main'
include { MINIMAP2_INDEX }          from '../modules/nf-core/minimap2/index/main'
include { MINIMAP2_ALIGN }          from '../modules/nf-core/minimap2/align/main'
include { SAMTOOLS_SORT }           from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_COVERAGE }       from '../modules/nf-core/samtools/coverage/main'
include { SEMIBIN_SINGLEEASYBIN }   from '../modules/nf-core/semibin/singleeasybin/main'
include { CHECKM2_PREDICT }         from '../modules/nf-core/checkm2/predict/main'
include { CHECKM2_PARSE }           from '../modules/local/checkm2/parse/main'  // Use local fixed version
include { BAKTA_BAKTA }             from '../modules/local/bakta/bakta/main'
include { SINGLEM_PIPE }            from '../modules/local/singlem/pipe/main'
include { SINGLEM_PIPE as SINGLEM_PIPE_BINS } from '../modules/local/singlem/pipe/main'
include { SINGLEM_APPRAISE }        from '../modules/local/singlem/appraise/main'
include { TAXBURST }                from '../modules/local/taxburst/main'

// Define default parameters
params.input_dir = ''
params.checkm2_db = ''
params.bakta_db = ''
params.singlem_metapackage = '' // Path to SingleM metapackage
params.output_dir = './results'
params.flye_mode = '--nano-hq'
params.semibin_environment = 'human_gut'
params.completeness_threshold = 80.0  // Threshold for generating Bakta plots

workflow SOMATEM_MAGS {

    take:
    reads               // channel: [ val(meta), path(reads) ]
    checkm2_db          // channel: [ val(dbmeta), path(db) ]
    bakta_db            // channel: path(db)
    singlem_metapackage // channel: path(metapackage)
    flye_mode           // val: sequencing data type for Flye
    semibin_environment // val: sample environment for SemiBin2

    main:
    ch_versions = Channel.empty()

    // STEP 1: Taxonomic profiling with SingleM on raw reads (metagenome) - RUN FIRST
    ch_metagenome_reads = reads.map { meta, reads_file ->
        def new_meta = meta.clone()
        new_meta.single_end = true
        new_meta.input_type = 'reads'
        [new_meta, reads_file]
    }
    
    SINGLEM_PIPE(ch_metagenome_reads, singlem_metapackage)
    ch_versions = ch_versions.mix(SINGLEM_PIPE.out.versions)

    // STEP 2: Create interactive taxonomic visualization with TaxBurst - RUN SECOND
    TAXBURST(SINGLEM_PIPE.out.taxonomic_profile, 'SingleM')
    ch_versions = ch_versions.mix(TAXBURST.out.versions)

    // STEP 3: Assembly with Flye - RUN AFTER TAXONOMIC PROFILING
    FLYE(reads, flye_mode)
    ch_versions = ch_versions.mix(FLYE.out.versions)

    // Create minimap2 index
    MINIMAP2_INDEX(FLYE.out.fasta)
    ch_versions = ch_versions.mix(MINIMAP2_INDEX.out.versions)

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
    SAMTOOLS_SORT(MINIMAP2_ALIGN.out.bam, FLYE.out.fasta)
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions)

    // Calculate coverage
    ch_bam_for_coverage = SAMTOOLS_SORT.out.bam.map { meta, bam ->
        [meta, bam, []] // Add empty index slot
    }
    
    ch_fasta_for_coverage = FLYE.out.fasta
    
    ch_fqi_for_coverage = SAMTOOLS_SORT.out.bam.map { meta, bam ->
        [meta, file('OPTIONAL_FILE')]
    }
    
    SAMTOOLS_COVERAGE(
        ch_bam_for_coverage,
        ch_fasta_for_coverage,
        ch_fqi_for_coverage
    )
    ch_versions = ch_versions.mix(SAMTOOLS_COVERAGE.out.versions)

    // Binning with SemiBin2
    ch_asm_bam = FLYE.out.fasta.join(SAMTOOLS_SORT.out.bam, by: [0])
    
    // Set SemiBin2 environment parameter
    ch_asm_bam_with_env = ch_asm_bam.map { meta, fasta, bam ->
        def new_meta = meta.clone()
        new_meta.semibin_env = semibin_environment
        [new_meta, fasta, bam]
    }
    
    SEMIBIN_SINGLEEASYBIN(ch_asm_bam_with_env)
    ch_versions = ch_versions.mix(SEMIBIN_SINGLEEASYBIN.out.versions)

    // Quality assessment with CheckM2
    CHECKM2_PREDICT(SEMIBIN_SINGLEEASYBIN.out.output_fasta, checkm2_db)
    ch_versions = ch_versions.mix(CHECKM2_PREDICT.out.versions)

    // Parse CheckM2 results to get completeness information
    CHECKM2_PARSE(CHECKM2_PREDICT.out.checkm2_tsv, SEMIBIN_SINGLEEASYBIN.out.output_fasta)
    ch_versions = ch_versions.mix(CHECKM2_PARSE.out.versions)

    // Run SingleM pipe on bins
    ch_bins_for_singlem = SEMIBIN_SINGLEEASYBIN.out.output_fasta.map { meta, bins ->
        def new_meta = meta.clone()
        new_meta.input_type = 'genome'
        [new_meta, bins]
    }
    
    SINGLEM_PIPE_BINS(ch_bins_for_singlem, singlem_metapackage)
    ch_versions = ch_versions.mix(SINGLEM_PIPE_BINS.out.versions)

    // SingleM appraise
    ch_metagenome_otu = SINGLEM_PIPE.out.otu_table
        .map { meta, otu -> 
            def clean_meta = [id: meta.id]
            return [clean_meta, otu]
        }
    
    ch_bins_otu = SINGLEM_PIPE_BINS.out.otu_table
        .map { meta, otu -> 
            def clean_meta = [id: meta.id]
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
            
            // Create staged files with unique names to avoid filename collision
            def appraise_dir = file("${params.output_dir}/appraise_staging")
            appraise_dir.mkdirs()
            
            def staged_metOtu = file("${appraise_dir}/${meta.id}_metagenome_otu_table.csv")
            def staged_binOtu = file("${appraise_dir}/${meta.id}_bins_otu_table.csv")
            
            // Copy files with unique names
            actualMetOtu.copyTo(staged_metOtu)
            actualBinOtu.copyTo(staged_binOtu)
            
            [new_meta, staged_metOtu, staged_binOtu, []]
        }

    SINGLEM_APPRAISE(ch_appraise_input, singlem_metapackage)
    ch_versions = ch_versions.mix(SINGLEM_APPRAISE.out.versions)

    // Enhanced Bakta annotation with completeness-based plot generation
    ch_bins_for_annotation = SEMIBIN_SINGLEEASYBIN.out.output_fasta
        .transpose()
        .map { meta, bin ->
            def new_meta = meta.clone()
            new_meta.id = "${meta.id}_${bin.baseName}"
            
            // Extract the SemiBin part from the filename for completeness lookup
            // e.g., s4_final_SemiBin_11.fa.gz -> SemiBin_11 or SemiBin_11.fa.gz -> SemiBin_11
            def bin_name = bin.name  // Get full filename with extension
            def semibin_name = bin_name.replaceAll(/.*?(SemiBin_\d+)\.fa(\.gz)?$/, '$1')
            new_meta.bin_name = semibin_name
            [new_meta, bin]
        }

    // Add completeness information to bin metadata
    ch_completeness_map = CHECKM2_PARSE.out.completeness_map
        .map { meta, csv_file ->
            // Parse the CSV file to create a map
            def completeness_map = [:]
            csv_file.text.split('\n').drop(1).each { line ->
                if (line.trim()) {
                    def parts = line.split(',')
                    if (parts.size() >= 2) {
                        completeness_map[parts[0]] = Float.parseFloat(parts[1])
                    }
                }
            }
            return completeness_map
        }

    // Combine bins with completeness information and filter for high-quality bins only
    ch_bins_with_completeness = ch_bins_for_annotation
        .combine(ch_completeness_map)
        .map { meta, bin, completeness_map ->
            def new_meta = meta.clone()
            def bin_name = meta.bin_name
            
            // Look up completeness for this bin
            def completeness = completeness_map[bin_name]
            if (completeness != null) {
                new_meta.completeness = completeness
                println "Bin ${bin_name}: ${completeness}% complete - ${completeness >= params.completeness_threshold ? 'WILL' : 'WILL NOT'} run Bakta annotation"
            } else {
                new_meta.completeness = null
                println "Bin ${bin_name}: completeness unknown - will not run Bakta annotation"
            }
            
            [new_meta, bin]
        }
        .filter { meta, bin ->
            // Only process bins that meet the completeness threshold
            return meta.completeness != null && meta.completeness >= params.completeness_threshold
        }

    // Only run Bakta on high-quality bins (≥completeness_threshold)
    BAKTA_BAKTA(ch_bins_with_completeness, bakta_db, [], [])
    ch_versions = ch_versions.mix(BAKTA_BAKTA.out.versions)

    // Note: File publishing is handled in the main workflow block to avoid subscribe issues

    emit:
    // Original outputs
    assembly        = FLYE.out.fasta
    assembly_gfa    = FLYE.out.gfa
    assembly_log    = FLYE.out.log
    bam_sorted      = SAMTOOLS_SORT.out.bam
    coverage        = SAMTOOLS_COVERAGE.out.coverage
    bins            = SEMIBIN_SINGLEEASYBIN.out.output_fasta
    bins_csv        = SEMIBIN_SINGLEEASYBIN.out.csv
    bins_tsv        = SEMIBIN_SINGLEEASYBIN.out.tsv
    checkm2_report  = CHECKM2_PREDICT.out.checkm2_tsv
    checkm2_output  = CHECKM2_PREDICT.out.checkm2_output
    completeness_map = CHECKM2_PARSE.out.completeness_map
    
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
    appraise_binned_otu = SINGLEM_APPRAISE.out.binned_otu_table
    appraise_unbinned_otu = SINGLEM_APPRAISE.out.unbinned_otu_table
    appraise_assembled_otu = SINGLEM_APPRAISE.out.assembled_otu_table
    appraise_unaccounted_otu = SINGLEM_APPRAISE.out.unaccounted_otu_table
    appraise_plot = SINGLEM_APPRAISE.out.plot
    
    versions        = ch_versions
}

// Main workflow for direct execution
workflow {
    // Validate required parameters
    if (!params.input_dir) {
        error "Please provide --input_dir parameter"
    }
    if (!params.checkm2_db) {
        error "Please provide --checkm2_db parameter"
    }
    if (!params.bakta_db) {
        error "Please provide --bakta_db parameter"
    }
    if (!params.singlem_metapackage) {
        error "Please provide --singlem_metapackage parameter"
    }

    // Create output directory
    file(params.output_dir).mkdirs()

    // Prepare input channels
    ch_reads = Channel.fromPath("${params.input_dir}/*.fastq.gz")
        .map { file -> 
            def meta = [id: file.baseName.replaceAll(/\.fastq(\.gz)?$/, '')]
            [meta, file]
        }
    
    ch_checkm2_db = Channel.fromPath(params.checkm2_db)
        .map { file -> 
            def meta = [id: file.name]
            [meta, file]
        }
    
    ch_bakta_db = Channel.value(params.bakta_db)
    ch_singlem_metapackage = Channel.value(params.singlem_metapackage)

    // Run the subworkflow
    SOMATEM_MAGS(
        ch_reads, 
        ch_checkm2_db, 
        ch_bakta_db, 
        ch_singlem_metapackage,
        params.flye_mode, 
        params.semibin_environment
    )

    // Display progress with completeness information
    SOMATEM_MAGS.out.singlem_profile.view { meta, profile -> "✓ Taxonomic profiling completed for ${meta.id}" }
    SOMATEM_MAGS.out.taxburst_html.view { meta, html -> "✓ Interactive taxonomy visualization created for ${meta.id}" }
    SOMATEM_MAGS.out.assembly.view { meta, fasta -> "✓ Assembly completed for ${meta.id}" }
    SOMATEM_MAGS.out.coverage.view { meta, coverage -> "✓ Coverage calculated for ${meta.id}" }
    SOMATEM_MAGS.out.bins.view { meta, bins -> "✓ Binning completed for ${meta.id}: ${bins.size()} bins" }
    SOMATEM_MAGS.out.checkm2_report.view { meta, report -> "✓ Quality assessment completed for ${meta.id}" }
    
    // Show which bins will get plots
    SOMATEM_MAGS.out.bakta_embl.view { meta, embl -> 
        def completeness = meta.completeness ?: "unknown"
        def plot_status = (meta.completeness != null && meta.completeness >= params.completeness_threshold) ? "with plots" : "without plots"
        "✓ Bakta annotation completed for ${meta.id} (${completeness}% complete) ${plot_status}"
    }
    
    SOMATEM_MAGS.out.appraise_summary.view { meta, summary -> "✓ SingleM appraise analysis completed for ${meta.id}" }

    // Count high-quality bins
    SOMATEM_MAGS.out.bakta_embl
        .map { meta, embl -> meta.completeness }
        .filter { it != null && it >= params.completeness_threshold }
        .count()
        .view { count -> "✓ Generated plots for ${count} high-quality bins (≥${params.completeness_threshold}% complete)" }

    // Publish results to organized directories with safe copying
    SOMATEM_MAGS.out.assembly.subscribe { meta, fasta ->
        def dest = file("${params.output_dir}/assembly/${meta.id}_assembly.fasta")
        dest.parent.mkdirs()
        try {
            if (fasta instanceof Collection && fasta.size() > 0) {
                fasta[0].copyTo(dest)
            } else {
                fasta.copyTo(dest)
            }
        } catch (Exception e) {
            println "Warning: Could not copy assembly for ${meta.id}: ${e.message}"
        }
    }

    SOMATEM_MAGS.out.assembly_gfa.subscribe { meta, gfa ->
        def dest = file("${params.output_dir}/assembly/${meta.id}_assembly.gfa")
        dest.parent.mkdirs()
        try {
            if (gfa instanceof Collection && gfa.size() > 0) {
                gfa[0].copyTo(dest)
            } else {
                gfa.copyTo(dest)
            }
        } catch (Exception e) {
            println "Warning: Could not copy assembly GFA for ${meta.id}: ${e.message}"
        }
    }

    SOMATEM_MAGS.out.assembly_log.subscribe { meta, log ->
        def dest = file("${params.output_dir}/assembly/${meta.id}_flye.log")
        dest.parent.mkdirs()
        try {
            if (log instanceof Collection && log.size() > 0) {
                log[0].copyTo(dest)
            } else {
                log.copyTo(dest)
            }
        } catch (Exception e) {
            println "Warning: Could not copy assembly log for ${meta.id}: ${e.message}"
        }
    }

    SOMATEM_MAGS.out.bam_sorted.subscribe { meta, bam ->
        def dest = file("${params.output_dir}/mapping/${meta.id}_sorted.bam")
        dest.parent.mkdirs()
        try {
            if (bam instanceof Collection && bam.size() > 0) {
                bam[0].copyTo(dest)
            } else {
                bam.copyTo(dest)
            }
        } catch (Exception e) {
            println "Warning: Could not copy BAM for ${meta.id}: ${e.message}"
        }
    }

    SOMATEM_MAGS.out.coverage.subscribe { meta, coverage ->
        def dest = file("${params.output_dir}/mapping/${meta.id}_coverage.txt")
        dest.parent.mkdirs()
        try {
            if (coverage instanceof Collection && coverage.size() > 0) {
                coverage[0].copyTo(dest)
            } else {
                coverage.copyTo(dest)
            }
        } catch (Exception e) {
            println "Warning: Could not copy coverage for ${meta.id}: ${e.message}"
        }
    }

    SOMATEM_MAGS.out.bins_csv.subscribe { meta, csv ->
        def dest = file("${params.output_dir}/binning/${meta.id}_semibin_results.csv")
        dest.parent.mkdirs()
        try {
            if (csv instanceof Collection && csv.size() > 0) {
                csv[0].copyTo(dest)
            } else {
                csv.copyTo(dest)
            }
        } catch (Exception e) {
            println "Warning: Could not copy bins CSV for ${meta.id}: ${e.message}"
        }
    }

    SOMATEM_MAGS.out.bins_tsv.subscribe { meta, tsv ->
        def dest = file("${params.output_dir}/binning/${meta.id}_semibin_results.tsv")
        dest.parent.mkdirs()
        try {
            if (tsv instanceof Collection && tsv.size() > 0) {
                tsv[0].copyTo(dest)
            } else {
                tsv.copyTo(dest)
            }
        } catch (Exception e) {
            println "Warning: Could not copy bins TSV for ${meta.id}: ${e.message}"
        }
    }

    SOMATEM_MAGS.out.checkm2_report.subscribe { meta, tsv ->
        def dest = file("${params.output_dir}/quality/${meta.id}_checkm2_report.tsv")
        dest.parent.mkdirs()
        try {
            if (tsv instanceof Collection && tsv.size() > 0) {
                tsv[0].copyTo(dest)
            } else {
                tsv.copyTo(dest)
            }
        } catch (Exception e) {
            println "Warning: Could not copy CheckM2 report for ${meta.id}: ${e.message}"
        }
    }

    SOMATEM_MAGS.out.completeness_map.subscribe { meta, csv ->
        def dest = file("${params.output_dir}/quality/${meta.id}_completeness_map.csv")
        dest.parent.mkdirs()
        try {
            if (csv instanceof Collection && csv.size() > 0) {
                csv[0].copyTo(dest)
            } else {
                csv.copyTo(dest)
            }
        } catch (Exception e) {
            println "Warning: Could not copy completeness map for ${meta.id}: ${e.message}"
        }
    }

    SOMATEM_MAGS.out.singlem_profile.subscribe { meta, profile ->
        def dest = file("${params.output_dir}/taxonomy/${meta.id}_metagenome_taxonomic_profile.csv")
        dest.parent.mkdirs()
        try {
            if (profile instanceof Collection && profile.size() > 0) {
                profile[0].copyTo(dest)
            } else {
                profile.copyTo(dest)
            }
        } catch (Exception e) {
            println "Warning: Could not copy SingleM profile for ${meta.id}: ${e.message}"
            println "Profile object type: ${profile.getClass()}"
        }
    }

    SOMATEM_MAGS.out.singlem_otu.subscribe { meta, otu ->
        def dest = file("${params.output_dir}/taxonomy/${meta.id}_metagenome_otu_table.csv")
        dest.parent.mkdirs()
        try {
            if (otu instanceof Collection && otu.size() > 0) {
                otu[0].copyTo(dest)
            } else {
                otu.copyTo(dest)
            }
        } catch (Exception e) {
            println "Warning: Could not copy SingleM OTU table for ${meta.id}: ${e.message}"
            println "OTU object type: ${otu.getClass()}"
        }
    }

    SOMATEM_MAGS.out.taxburst_html.subscribe { meta, html ->
        def dest = file("${params.output_dir}/taxonomy/${meta.id}_taxburst_visualization.html")
        dest.parent.mkdirs()
        try {
            if (html instanceof Collection && html.size() > 0) {
                html[0].copyTo(dest)
            } else {
                html.copyTo(dest)
            }
        } catch (Exception e) {
            println "Warning: Could not copy TaxBurst HTML for ${meta.id}: ${e.message}"
            println "HTML object type: ${html.getClass()}"
        }
    }

    // Handle bins output carefully since it's a collection
    SOMATEM_MAGS.out.bins.subscribe { meta, bins ->
        def dest_dir = file("${params.output_dir}/binning/${meta.id}_bins")
        dest_dir.mkdirs()
        try {
            if (bins instanceof Collection || bins instanceof List) {
                bins.each { bin ->
                    if (bin != null && bin.exists()) {
                        def dest_file = file("${dest_dir}/${bin.name}")
                        bin.copyTo(dest_file)
                    }
                }
            } else if (bins != null && bins.exists()) {
                def dest_file = file("${dest_dir}/${bins.name}")
                bins.copyTo(dest_file)
            }
        } catch (Exception e) {
            println "Warning: Could not copy bins for ${meta.id}: ${e.message}"
            println "Bins object type: ${bins.getClass()}"
            println "Bins content: ${bins}"
        }
    }

    // Publish all Bakta annotation outputs organized by bin
    // Helper function to safely copy files
    def copyBaktaFile = { meta, file_obj, extension ->
        if (file_obj && file_obj.exists()) {
            def bin_dir = file("${params.output_dir}/annotation/${meta.id}")
            bin_dir.mkdirs()
            def dest = file("${bin_dir}/${meta.id}.${extension}")
            try {
                if (file_obj instanceof Collection && file_obj.size() > 0) {
                    file_obj[0].copyTo(dest)
                } else {
                    file_obj.copyTo(dest)
                }
                return dest
            } catch (Exception e) {
                println "Warning: Could not copy Bakta ${extension.toUpperCase()} for ${meta.id}: ${e.message}"
                return null
            }
        }
        return null
    }

    // Publish all Bakta outputs
    SOMATEM_MAGS.out.bakta_embl.subscribe { meta, embl ->
        copyBaktaFile(meta, embl, 'embl')
    }

    SOMATEM_MAGS.out.bakta_faa.subscribe { meta, faa ->
        copyBaktaFile(meta, faa, 'faa')
    }

    SOMATEM_MAGS.out.bakta_ffn.subscribe { meta, ffn ->
        copyBaktaFile(meta, ffn, 'ffn')
    }

    SOMATEM_MAGS.out.bakta_fna.subscribe { meta, fna ->
        copyBaktaFile(meta, fna, 'fna')
    }

    SOMATEM_MAGS.out.bakta_gbff.subscribe { meta, gbff ->
        copyBaktaFile(meta, gbff, 'gbff')
    }

    SOMATEM_MAGS.out.bakta_gff.subscribe { meta, gff ->
        copyBaktaFile(meta, gff, 'gff')
    }

    SOMATEM_MAGS.out.bakta_hypotheticals_tsv.subscribe { meta, hyp_tsv ->
        copyBaktaFile(meta, hyp_tsv, 'hypotheticals.tsv')
    }

    SOMATEM_MAGS.out.bakta_hypotheticals_faa.subscribe { meta, hyp_faa ->
        copyBaktaFile(meta, hyp_faa, 'hypotheticals.faa')
    }

    SOMATEM_MAGS.out.bakta_tsv.subscribe { meta, tsv ->
        copyBaktaFile(meta, tsv, 'tsv')
    }

    SOMATEM_MAGS.out.bakta_txt.subscribe { meta, txt ->
        copyBaktaFile(meta, txt, 'txt')
    }

    SOMATEM_MAGS.out.bakta_inference_tsv.subscribe { meta, inf_tsv ->
        copyBaktaFile(meta, inf_tsv, 'inference.tsv')
    }

    SOMATEM_MAGS.out.bakta_png.subscribe { meta, png ->
        copyBaktaFile(meta, png, 'png')
    }

    SOMATEM_MAGS.out.bakta_svg.subscribe { meta, svg ->
        copyBaktaFile(meta, svg, 'svg')
    }

    SOMATEM_MAGS.out.bakta_json.subscribe { meta, json ->
        copyBaktaFile(meta, json, 'json')
    }

    // Publish SingleM appraise outputs
    SOMATEM_MAGS.out.appraise_summary.subscribe { meta, summary ->
        def dest = file("${params.output_dir}/appraise/${meta.id}_summary.tsv")
        dest.parent.mkdirs()
        try {
            if (summary instanceof Collection && summary.size() > 0) {
                summary[0].copyTo(dest)
            } else {
                summary.copyTo(dest)
            }
        } catch (Exception e) {
            println "Warning: Could not copy SingleM appraise summary for ${meta.id}: ${e.message}"
        }
    }

    SOMATEM_MAGS.out.appraise_plot.subscribe { meta, plot ->
        if (plot && plot.exists()) {
            def dest = file("${params.output_dir}/appraise/${meta.id}_plot.svg")
            dest.parent.mkdirs()
            try {
                if (plot instanceof Collection && plot.size() > 0) {
                    plot[0].copyTo(dest)
                } else {
                    plot.copyTo(dest)
                }
                println "✓ Saved SingleM appraise plot: ${dest}"
            } catch (Exception e) {
                println "Warning: Could not copy SingleM appraise plot for ${meta.id}: ${e.message}"
            }
        }
    }
}