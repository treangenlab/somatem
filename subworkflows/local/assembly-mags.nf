#!/usr/bin/env nextflow

// Assembly_MAGS subworkflow: de novo assembly, mapping, binning, quality assessment, and annotation

// Main workflow for direct execution
workflow {
    
    // Prepare input channels
    ch_reads = Channel.fromPath("${params.input_dir}/*.fastq.gz")
        .ifEmpty { error "No .fastq.gz files found in ${params.input_dir}" }
        .map { file -> 
            def meta = [id: file.baseName.replaceAll(/\.fastq(\.gz)?$/, '')]
            println "Processing input file: ${file} -> sample ID: ${meta.id}"
            [meta, file]
        }
    
    // Create value channels for databases
    ch_checkm2_db = Channel.value(params.checkm2_db)
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

    // Display progress information
    SOMATEM_MAGS.out.singlem_profile.view { meta, profile -> "✓ Taxonomic profiling completed for ${meta.id}" }
    SOMATEM_MAGS.out.taxburst_html.view { meta, html -> "✓ Interactive taxonomy visualization created for ${meta.id}" }
    SOMATEM_MAGS.out.assembly.view { meta, fasta -> "✓ Assembly completed for ${meta.id}" }
    SOMATEM_MAGS.out.coverage.view { meta, coverage -> "✓ Coverage calculated for ${meta.id}" }
    SOMATEM_MAGS.out.bins.view { meta, bins -> "✓ Binning completed for ${meta.id}: ${bins.size()} bins" }
    SOMATEM_MAGS.out.checkm2_report.view { meta, report -> "✓ Quality assessment completed for ${meta.id}" }
    
    // Show which bins got annotation
    SOMATEM_MAGS.out.bakta_embl.view { meta, embl -> 
        def completeness = meta.completeness ?: "unknown"
        "✓ Bakta annotation completed for ${meta.id} (${completeness}% complete)"
    }
    
    SOMATEM_MAGS.out.appraise_summary.view { meta, summary -> "✓ SingleM appraise analysis completed for ${meta.id}" }

    // Count high-quality bins
    SOMATEM_MAGS.out.bakta_embl
        .map { meta, embl -> meta.completeness }
        .filter { it != null && it >= params.completeness_threshold }
        .count()
        .view { count -> "✓ Generated annotations for ${count} high-quality bins (≥${params.completeness_threshold}% complete)" }
}