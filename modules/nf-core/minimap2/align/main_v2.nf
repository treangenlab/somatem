/*
 * MINIMAP2_ALIGN MODULE
 * Aligns long reads to reference genome using minimap2
 * Based on marlin.py map_reads_minimap2() function
 */

process MINIMAP2_ALIGN {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/env.yml"

    input:
    tuple val(meta), path(reads), path(reference)

    output:
    tuple val(meta), path("*.bam"), path("*.bam.bai"), path(reference), emit: bam
    path "versions.yml", emit: versions

    script:
    def prefix = "${meta.id}"
    def samtools_threads = Math.max(1, (task.cpus / 2) as int)
    """
    # Map reads to reference with minimap2 (map-ont preset)
    minimap2 \\
        -x map-ont \\
        -a \\
        -t ${task.cpus} \\
        ${reference} \\
        ${reads} \\
        | samtools sort -@ ${samtools_threads} -o ${prefix}.bam -

    # Index the BAM file
    samtools index ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = "${meta.id}"
    """
    touch ${prefix}.bam
    touch ${prefix}.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: stub
        samtools: stub
    END_VERSIONS
    """
}
