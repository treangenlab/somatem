process SYLPH_PROFILE {
    tag "${meta.id}"
    label 'process_high'

    publishDir { "${params.outdir}/taxonomic_profiling/${meta.id}/sylph" },
        mode: params.publish_dir_mode,
        pattern: '*.tsv'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/sylph:0.9.0--ha6fb395_0'
        : 'biocontainers/sylph:0.9.0--ha6fb395_0'}"

    input:
    tuple val(meta), path(reads)
    path database

    output:
    tuple val(meta), path('*.tsv'), emit: profile_out
    path 'versions.yml', emit: versions
    tuple val("${task.process}"), val('sylph'), eval('sylph -V | awk "{print \$2}"'), topic: versions, emit: versions_sylph

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input = meta.single_end ? "${reads}" : "-1 ${reads[0]} -2 ${reads[1]}"
    """
    sylph profile \\
        -t ${task.cpus} \\
        ${args} \\
        ${database}\\
        ${input} \\
        -o ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sylph: \$(sylph -V | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sylph: "stub"
    END_VERSIONS
    """
}
