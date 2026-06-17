process SEQSCREEN_FORMAT_REPORT {
    tag "$meta.id"
    label 'process_single'

    conda "${projectDir}/modules/local/seqscreen/environment.yml"

    publishDir { "${params.outdir}/seqscreen/${meta.id}" },
        mode: params.publish_dir_mode,
        pattern: "seqscreen_output"

    input:
    tuple val(meta), path(report), path(html), path(taxonomy), path(functional), path(reference_inference)
    path db
    path assets
    val mode

    output:
    tuple val(meta), path("seqscreen_output"), emit: outdir
    tuple val(meta), path("seqscreen_output/report_generation/*seqscreen_report.tsv"), emit: report
    tuple val(meta), path("seqscreen_output/report_generation/seqscreen_html_report"), emit: html
    tuple val(meta), path("seqscreen_output/taxonomic_identification/taxonomic_assignment/taxonomic_results.txt"), emit: taxonomic_results
    tuple val(meta), path("seqscreen_output/functional_annotation/functional_assignments/functional_results.txt"), emit: functional_results
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def format = params.seqscreen_format ?: 1
    def confidence = params.seqscreen_taxonomy_confidence_threshold ?: 0
    def filter_taxon = params.seqscreen_filter_taxon ?: '""'
    def keep_taxon = params.seqscreen_keep_taxon ?: '""'
    """
    mkdir -p seqscreen_output/report_generation
    mkdir -p seqscreen_output/taxonomic_identification/taxonomic_assignment
    mkdir -p seqscreen_output/functional_annotation/functional_assignments
    mkdir -p seqscreen_output/taxonomic_identification/taxonomic_assignment/inference_working

    cp ${report} seqscreen_output/report_generation/
    cp -R ${html} seqscreen_output/report_generation/seqscreen_html_report
    cp ${taxonomy} seqscreen_output/taxonomic_identification/taxonomic_assignment/taxonomic_results.txt
    cp ${functional} seqscreen_output/functional_annotation/functional_assignments/functional_results.txt
    cp -R ${reference_inference} seqscreen_output/reference_inference

    python3 ${assets}/scripts/format.py \\
        --report seqscreen_output/report_generation \\
        --format ${format} \\
        --mode ${mode} \\
        --databases=${db} \\
        --taxonomy_confidence_threshold=${confidence} \\
        --filter-taxon ${filter_taxon} \\
        --keep-taxon ${keep_taxon}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqscreen: \$(seqscreen --version 2>&1 | sed 's/^SeqScreen v//; s/^.* //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p seqscreen_output/report_generation/seqscreen_html_report
    mkdir -p seqscreen_output/taxonomic_identification/taxonomic_assignment
    mkdir -p seqscreen_output/functional_annotation/functional_assignments
    touch seqscreen_output/report_generation/seqscreen_report.tsv
    touch seqscreen_output/report_generation/seqscreen_html_report/index.html
    touch seqscreen_output/taxonomic_identification/taxonomic_assignment/taxonomic_results.txt
    touch seqscreen_output/functional_annotation/functional_assignments/functional_results.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqscreen: "stub"
    END_VERSIONS
    """
}
