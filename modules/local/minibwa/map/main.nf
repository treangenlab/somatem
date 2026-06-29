process MINIBWA_MAP {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(query)
    path index

    output:
    tuple val(meta), path("*.paf"), emit: paf
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    minibwa map \\
        -f \\
        -b cs \\
        -t ${task.cpus} \\
        ${args} \\
        ${index}/reference.fasta \\
        ${query} \\
        > ${prefix}.paf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minibwa: \$(minibwa 2>&1 | head -n 1 | sed 's/^.*minibwa/minibwa/')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cat <<-END_PAF > ${prefix}.paf
    ${meta.id}	32	0	32	+	ref	32	0	32	31	32	60	cs:Z::20*at:11
    END_PAF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minibwa: "stub"
    END_VERSIONS
    """
}
