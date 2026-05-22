process POLYPOLISH_POLISH {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"

    publishDir { "${params.outdir}/isolate_analysis/${meta.id}/polishing/polypolish" },
        mode: params.publish_dir_mode,
        pattern: "*.polypolish.fasta"

    input:
    tuple val(meta), path(assembly), path(filtered_1), path(filtered_2)

    output:
    tuple val(meta), path("*.polypolish.fasta"), emit: assembly
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = params.polypolish_args ?: ''
    def prefix = meta.id
    """
    polypolish polish \\
        ${args} \\
        ${assembly} \\
        ${filtered_1} \\
        ${filtered_2} \\
        > ${prefix}.polypolish.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        polypolish: \$(polypolish --version 2>&1 | sed 's/^polypolish //')
    END_VERSIONS
    """
}
