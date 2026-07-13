process SEQSCREEN_MUMMER {
    tag "$meta.id"
    label 'process_medium'

    conda "${projectDir}/modules/local/seqscreen/mummer/environment.yml"

    input:
    tuple val(meta), path(fasta)
    path db

    output:
    tuple val(meta), path("*.mummer_re.txt"), emit: results
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = "${meta.id}"
    """
    mummer \\
        -maxmatch \\
        -l 4 \\
        ${fasta} \\
        ${db}/rebase/rebase.fna \\
        > ${prefix}.mummer_re.txt \\
        2> mummer.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mummer: "\$(nucmer --version 2>&1 | head -n 1 | sed 's/^.* //')"
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.mummer_re.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mummer: "stub"
    END_VERSIONS
    """
}
