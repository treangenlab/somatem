process COVERM {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::coverm=0.7.0"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(contigs)

    output:
    tuple val(meta), path("*.bam")        , emit: bam
    tuple val(meta), path("*.coverage.tsv"), emit: coverage
    tuple val(meta), path("*.summary.tsv") , emit: summary
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    coverm contig \\
        --reference ${contigs} \\
        --reads ${reads} \\
        --output-file ${prefix}.coverage.tsv \\
        --bam-file-cache-directory . \\
        --threads ${task.cpus} \\
        ${args}

    # Generate summary statistics
    coverm contig \\
        --reference ${contigs} \\
        --reads ${reads} \\
        --methods mean \\
        --output-file ${prefix}.summary.tsv \\
        --threads ${task.cpus}

    # Move BAM files to expected output
    mv *.bam ${prefix}.bam || true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coverm: \$(coverm --version 2>&1 | sed 's/CoverM //g')
    END_VERSIONS
    """
}