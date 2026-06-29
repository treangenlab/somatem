process MINIBWA_INDEX {
    tag "$reference"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"

    input:
    path reference

    output:
    path "minibwa_index", emit: index
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mkdir -p minibwa_index
    cp ${reference} minibwa_index/reference.fasta
    minibwa index \\
        -t ${task.cpus} \\
        ${args} \\
        minibwa_index/reference.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minibwa: \$(minibwa 2>&1 | head -n 1 | sed 's/^.*minibwa/minibwa/')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p minibwa_index
    cp ${reference} minibwa_index/reference.fasta
    touch minibwa_index/reference.fasta.l2b
    touch minibwa_index/reference.fasta.mbw

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minibwa: "stub"
    END_VERSIONS
    """
}
