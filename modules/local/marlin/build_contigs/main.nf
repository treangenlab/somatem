/*
 * BUILD_CONTIGS MODULE
 * Breaks consensus sequences at coverage gaps to produce contigs
 * Uses the build_contigs.py script from bin/
 */

process BUILD_CONTIGS {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(consensus), path(depth_file)
    val min_cov
    val min_gap_len
    val min_contig_len

    output:
    tuple val(meta), path("*.contigs.fasta"), emit: contigs
    tuple val(meta), path("*.layout.tsv"), emit: layout
    path "versions.yml", emit: versions

    script:
    def prefix = "${meta.id}"
    """
    python $PWD/bin/build_contigs.py \\
        --consensus ${consensus} \\
        --depth ${depth_file} \\
        --min-cov ${min_cov} \\
        --min-gap-len ${min_gap_len} \\
        --min-contig-len ${min_contig_len} \\
        --contigs-out ${prefix}.contigs.fasta \\
        --layout-out ${prefix}.layout.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def prefix = "${meta.id}"
    """
    touch ${prefix}.contigs.fasta
    touch ${prefix}.layout.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: stub
    END_VERSIONS
    """
}
