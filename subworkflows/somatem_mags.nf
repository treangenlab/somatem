#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Metagenomics analysis subworkflow: assembly, mapping, binning, quality assessment, and annotation

// Include nf-core modules
include { FLYE }                    from '../modules/nf-core/flye/main'
include { MINIMAP2_INDEX }          from '../modules/nf-core/minimap2/index/main'
include { MINIMAP2_ALIGN }          from '../modules/nf-core/minimap2/align/main'
include { SAMTOOLS_SORT }           from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX }          from '../modules/nf-core/samtools/index/main'
include { SAMTOOLS_COVERAGE }       from '../modules/nf-core/samtools/coverage/main'
include { SEMIBIN_SINGLEEASYBIN }   from '../modules/nf-core/semibin/singleeasybin/main'
include { CHECKM2_PREDICT }         from '../modules/nf-core/checkm2/predict/main'
include { BAKTA_BAKTA }             from '../modules/nf-core/bakta/bakta/main'

// Define default parameters
params.input_dir = ''
params.checkm2_db = ''
params.bakta_db = ''
params.output_dir = './results'
params.flye_mode = '--nano-hq'
params.semibin_environment = 'human_gut'

workflow SOMATEM_MAGS {

    take:
    reads               // channel: [ val(meta), path(reads) ]
    checkm2_db          // channel: [ val(dbmeta), path(db) ]
    bakta_db            // channel: path(db)
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

    // Index sorted BAM files
    SAMTOOLS_INDEX(SAMTOOLS_SORT.out.bam)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

    // Calculate coverage - FIXED
    ch_bam_idx = SAMTOOLS_SORT.out.bam.join(SAMTOOLS_INDEX.out.bai, by: [0])
    
    // Debug channels to see what's happening
    ch_bam_idx.view { "BAM+IDX: ${it}" }
    FLYE.out.fasta.view { "FASTA: ${it}" }
    
    SAMTOOLS_COVERAGE(
        ch_bam_idx, 
        FLYE.out.fasta, 
        Channel.empty()
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
    assembly        = FLYE.out.fasta
    assembly_gfa    = FLYE.out.gfa
    assembly_log    = FLYE.out.log
    bam_sorted      = SAMTOOLS_SORT.out.bam
    bam_index       = SAMTOOLS_INDEX.out.bai
    coverage        = SAMTOOLS_COVERAGE.out.coverage
    bins            = SEMIBIN_SINGLEEASYBIN.out.output_fasta
    bins_csv        = SEMIBIN_SINGLEEASYBIN.out.csv
    bins_tsv        = SEMIBIN_SINGLEEASYBIN.out.tsv
    checkm2_report  = CHECKM2_PREDICT.out.checkm2_tsv
    checkm2_output  = CHECKM2_PREDICT.out.checkm2_output
    bakta_embl      = BAKTA_BAKTA.out.embl
    bakta_faa       = BAKTA_BAKTA.out.faa
    bakta_gbff      = BAKTA_BAKTA.out.gbff
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

    // Run the subworkflow
    SOMATEM_MAGS(
        ch_reads, 
        ch_checkm2_db, 
        ch_bakta_db, 
        params.flye_mode, 
        params.semibin_environment
    )

    // Publish results to output directory with explicit copying
    SOMATEM_MAGS.out.assembly
        .subscribe { meta, fasta ->
            def dest = file("${params.output_dir}/assembly/${meta.id}.fasta")
            dest.parent.mkdirs()
            fasta.copyTo(dest)
            println "Assembly saved: ${dest}"
        }

    SOMATEM_MAGS.out.coverage
        .subscribe { meta, coverage ->
            def dest = file("${params.output_dir}/coverage/${meta.id}_coverage.txt")
            dest.parent.mkdirs()
            coverage.copyTo(dest)
            println "Coverage saved: ${dest}"
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

    // Display progress
    SOMATEM_MAGS.out.assembly.view { meta, fasta -> "✓ Assembly completed for ${meta.id}" }
    SOMATEM_MAGS.out.coverage.view { meta, coverage -> "✓ Coverage calculated for ${meta.id}" }
    SOMATEM_MAGS.out.bins.view { meta, bins -> "✓ Binning completed for ${meta.id}: ${bins.size()} bins" }
    SOMATEM_MAGS.out.checkm2_report.view { meta, report -> "✓ Quality assessment completed for ${meta.id}" }
}