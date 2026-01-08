/*
 * SAMTOOLS_DEPTH MODULE
 * Calculates depth of coverage from BAM file
 * Based on marlin.py compute_depth() function
 */

process SAMTOOLS_DEPTH {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(bam), path(bai), path(reference)

    output:
    tuple val(meta), path("*.depth.txt"), emit: depth
    path "versions.yml", emit: versions

    script:
    def prefix = "${meta.id}"
    """
    # Compute per-base depth with -a flag (all positions including zero coverage)
    samtools depth -a ${bam} > ${prefix}.depth.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = "${meta.id}"
    """
    touch ${prefix}.depth.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: stub
    END_VERSIONS
    """
}
