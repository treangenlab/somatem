process SEQSCREEN_CENTRIFUGE {
    tag "$meta.id"
    label 'process_medium'

    conda "${projectDir}/modules/local/seqscreen/centrifuge/environment.yml"

    input:
    tuple val(meta), path(fasta)
    path db

    output:
    tuple val(meta), path("*.centrifuge"), emit: results
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    centrifuge \\
        -x ${db}/centrifuge/abv \\
        -f \\
        -U ${fasta} \\
        -S ${prefix}.centrifuge \\
        --threads ${task.cpus} \\
        --seed 1 \\
        > centrifuge.log 2>&1

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        centrifuge: "\$(centrifuge --version 2>&1 | head -n 1 | sed 's/^.*centrifuge-class version //')"
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.centrifuge
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        centrifuge: "stub"
    END_VERSIONS
    """
}
