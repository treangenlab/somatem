process SEQSCREEN_WINDOW_FASTA {
    tag "$meta.id"
    label 'process_single'

    conda "${projectDir}/modules/local/seqscreen/window_fasta/environment.yml"

    input:
    tuple val(meta), path(fasta)
    path assets

    output:
    tuple val(meta), path("*.windowed.fasta"), emit: fasta
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def length = params.seqscreen_window_length ?: 200
    def overlap = params.seqscreen_window_overlap ?: 100
    """
    python3 ${assets}/scripts/split_window_input.py \\
        --fasta ${fasta} \\
        --output ${prefix}.windowed.fasta \\
        --length ${length} \\
        --overlap ${overlap}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqscreen: "4.5"
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.windowed.fasta
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqscreen: "stub"
    END_VERSIONS
    """
}
