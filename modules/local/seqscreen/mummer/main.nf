process SEQSCREEN_MUMMER {
    tag "$meta.id"
    label 'process_medium'

    conda "${projectDir}/modules/local/seqscreen/environment.yml"

    input:
    tuple val(meta), path(fasta)
    path db
    path assets

    output:
    tuple val(meta), path("*.mummer_re.txt"), emit: results
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = "${meta.id}"
    """
    ${assets}/modules/mummer.sh \\
        --fasta=${fasta} \\
        --database=${db}/rebase/rebase.fna \\
        --out=${prefix}.mummer_re.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mummer: \$(nucmer --version 2>&1 | head -n 1 | sed 's/^.* //')
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
