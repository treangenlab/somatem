#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Metagenomics analysis subworkflow: assembly, mapping, binning, quality assessment, annotation, and taxonomic profiling

// Include nf-core modules
include { FLYE }                    from '../modules/nf-core/flye/main'
include { MINIMAP2_INDEX }          from '../modules/nf-core/minimap2/index/main'
include { MINIMAP2_ALIGN }          from '../modules/nf-core/minimap2/align/main'
include { SAMTOOLS_SORT }           from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_COVERAGE }       from '../modules/nf-core/samtools/coverage/main'
include { SEMIBIN_SINGLEEASYBIN }   from '../modules/nf-core/semibin/singleeasybin/main'
include { CHECKM2_PREDICT }         from '../modules/nf-core/checkm2/predict/main'
include { BAKTA_BAKTA }             from '../modules/nf-core/bakta/bakta/main'
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

    // Assembly with Flye
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

    // Calculate coverage - Handle the actual SAMTOOLS_SORT output structure
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

    // Taxonomic profiling with SingleM on raw reads (metagenome)
    ch_metagenome_reads = reads.map { meta, reads_file ->
        def new_meta = meta.clone()
        new_meta.single_end = true
        new_meta.input_type = 'reads'
        [new_meta, reads_file]
    }
    
    SINGLEM_PIPE(ch_metagenome_reads, singlem_metapackage)
    ch_versions = ch_versions.mix(SINGLEM_PIPE.out.versions)

    // Create interactive taxonomic visualization with TaxBurst
    TAXBURST(SINGLEM_PIPE.out.taxonomic_profile, 'SingleM')
    ch_versions = ch_versions.mix(TAXBURST.out.versions)

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

    // Run SingleM pipe on bins - Keep the SAME meta.id for successful joining
    ch_bins_for_singlem = SEMIBIN_SINGLEEASYBIN.out.output_fasta.map { meta, bins ->
        // Use the exact same meta structure as the metagenome channel
        def new_meta = meta.clone()  // Keep all original meta properties
        new_meta.input_type = 'genome'  // Just add the input_type
        [new_meta, bins]
    }
    
    // Run SingleM pipe on bins using the aliased module
    SINGLEM_PIPE_BINS(ch_bins_for_singlem, singlem_metapackage)
    ch_versions = ch_versions.mix(SINGLEM_PIPE_BINS.out.versions)

    // FIXED: SingleM appraise - Resolve filename collision using Nextflow's staging
    ch_metagenome_otu = SINGLEM_PIPE.out.otu_table
        .map { meta, otu -> 
            println "DEBUG: Metagenome OTU - meta.id: ${meta.id}, file: ${otu.name}"
            // Use Nextflow's staging to rename the file
            def clean_meta = [id: meta.id]
            return [clean_meta, otu]
        }
    
    ch_bins_otu = SINGLEM_PIPE_BINS.out.otu_table
        .map { meta, otu -> 
            println "DEBUG: Bins OTU - meta.id: ${meta.id}, file: ${otu.name}"
            // Use Nextflow's staging to rename the file
            def clean_meta = [id: meta.id]
            return [clean_meta, otu]
        }

    // Use cross product to ensure both channels are processed
    ch_appraise_input = ch_metagenome_otu
        .cross(ch_bins_otu) { it[0].id }  // Cross by meta.id
        .map { metagenome_tuple, bins_tuple ->
            def meta = metagenome_tuple[0]
            def metOtu = metagenome_tuple[1]
            def binOtu = bins_tuple[1]
            
            def new_meta = [id: "${meta.id}_appraisal"]
            println "DEBUG: Creating appraise input - meta: ${new_meta.id}"
            println "  - Metagenome OTU: ${metOtu.name}"
            println "  - Bins OTU: ${binOtu.name}"
            
            // FIXED: Create the renamed OTU files in the appraise output directory
            def appraise_dir = file("${params.output_dir}/appraise")
            appraise_dir.mkdirs()
            
            def staged_metOtu = file("${appraise_dir}/${meta.id}_metagenome_otu_table.csv")
            def staged_binOtu = file("${appraise_dir}/${meta.id}_bins_otu_table.csv")
            
            metOtu.copyTo(staged_metOtu)
            binOtu.copyTo(staged_binOtu)
            
            println "  - Staged Metagenome OTU: ${staged_metOtu}"
            println "  - Staged Bins OTU: ${staged_binOtu}"
            
            // SINGLEM_APPRAISE expects: metagenome_otu_tables, genome_otu_tables, assembly_otu_tables
            [new_meta, staged_metOtu, staged_binOtu, []]
        }
        .view { "DEBUG: Final appraise input: $it" }

    // Run SINGLEM_APPRAISE
    SINGLEM_APPRAISE(ch_appraise_input, singlem_metapackage)
    ch_versions = ch_versions.mix(SINGLEM_APPRAISE.out.versions)

    // Annotation with Bakta - transpose bins for individual annotation
    ch_bins_for_annotation = SEMIBIN_SINGLEEASYBIN.out.output_fasta
        .transpose()
        .map { meta, bin ->
            def new_meta = meta.clone()
            new_meta.id = "${meta.id}_${bin.baseName}"
            [new_meta, bin]
        }

    BAKTA_BAKTA(ch_bins_for_annotation, bakta_db, [], [])
    ch_versions = ch_versions.mix(BAKTA_BAKTA.out.versions)

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
    bakta_embl      = BAKTA_BAKTA.out.embl
    bakta_faa       = BAKTA_BAKTA.out.faa
    bakta_gbff      = BAKTA_BAKTA.out.gbff
    
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

    // Publish results with proper null checking
    SOMATEM_MAGS.out.assembly
        .subscribe { meta, fasta ->
            def dest = file("${params.output_dir}/assembly/${meta.id}.fasta")
            dest.parent.mkdirs()
            fasta.copyTo(dest)
            println "Assembly saved: ${dest}"
        }

    SOMATEM_MAGS.out.coverage
        .subscribe { meta, coverage ->
            if (coverage) {
                def dest = file("${params.output_dir}/coverage/${meta.id}_coverage.txt")
                dest.parent.mkdirs()
                coverage.copyTo(dest)
                println "Coverage saved: ${dest}"
            }
        }

    SOMATEM_MAGS.out.bins
        .subscribe { meta, bins ->
            def dest_dir = file("${params.output_dir}/binning/${meta.id}")
            dest_dir.mkdirs()
            bins.each { bin ->
                bin.copyTo("${dest_dir}/${bin.name}")
            }
            println "Bins saved: ${dest_dir} (${bins.size()} bins)"
        }

    SOMATEM_MAGS.out.checkm2_report
        .subscribe { meta, report ->
            def dest = file("${params.output_dir}/quality/${meta.id}_checkm2.tsv")
            dest.parent.mkdirs()
            report.copyTo(dest)
            println "CheckM2 report saved: ${dest}"
        }

    // Taxonomic profiling outputs
    SOMATEM_MAGS.out.singlem_profile
        .subscribe { meta, profile ->
            def dest = file("${params.output_dir}/taxonomy/${meta.id}_singlem.profile.tsv")
            dest.parent.mkdirs()
            profile.copyTo(dest)
            println "SingleM profile saved: ${dest}"
        }

    SOMATEM_MAGS.out.taxburst_html
        .subscribe { meta, html ->
            def dest = file("${params.output_dir}/taxonomy/${meta.id}_taxburst.html")
            dest.parent.mkdirs()
            html.copyTo(dest)
            println "TaxBurst visualization saved: ${dest}"
        }

    // SingleM appraise outputs - Fixed null handling
    SOMATEM_MAGS.out.appraise_summary
        .filter { meta, summary -> summary != null }
        .subscribe { meta, summary ->
            def dest = file("${params.output_dir}/appraise/${meta.id}_summary.txt")
            dest.parent.mkdirs()
            summary.copyTo(dest)
            println "SingleM appraise summary saved: ${dest}"
        }

    SOMATEM_MAGS.out.appraise_binned_otu
        .filter { meta, binned -> binned != null }
        .subscribe { meta, binned ->
            def dest = file("${params.output_dir}/appraise/${meta.id}_binned.csv")
            dest.parent.mkdirs()
            binned.copyTo(dest)
            println "SingleM appraise binned OTU table saved: ${dest}"
        }

    SOMATEM_MAGS.out.appraise_unbinned_otu
        .filter { meta, unbinned -> unbinned != null }
        .subscribe { meta, unbinned ->
            def dest = file("${params.output_dir}/appraise/${meta.id}_unbinned.csv")
            dest.parent.mkdirs()
            unbinned.copyTo(dest)
            println "SingleM appraise unbinned OTU table saved: ${dest}"
        }

    SOMATEM_MAGS.out.appraise_plot
        .filter { meta, plot -> plot != null }
        .subscribe { meta, plot ->
            def dest = file("${params.output_dir}/appraise/${meta.id}_plot.svg")
            dest.parent.mkdirs()
            plot.copyTo(dest)
            println "SingleM appraise plot saved: ${dest}"
        }

    // Display progress
    SOMATEM_MAGS.out.assembly.view { meta, fasta -> "✓ Assembly completed for ${meta.id}" }
    SOMATEM_MAGS.out.coverage.view { meta, coverage -> "✓ Coverage calculated for ${meta.id}" }
    SOMATEM_MAGS.out.bins.view { meta, bins -> "✓ Binning completed for ${meta.id}: ${bins.size()} bins" }
    SOMATEM_MAGS.out.checkm2_report.view { meta, report -> "✓ Quality assessment completed for ${meta.id}" }
    SOMATEM_MAGS.out.singlem_profile.view { meta, profile -> "✓ Taxonomic profiling completed for ${meta.id}" }
    SOMATEM_MAGS.out.taxburst_html.view { meta, html -> "✓ Interactive taxonomy visualization created for ${meta.id}" }
    SOMATEM_MAGS.out.appraise_summary.view { meta, summary -> "✓ SingleM appraise analysis completed for ${meta.id}" }
}