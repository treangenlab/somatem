process SEQSCREEN_VALIDATE_FASTA {
    tag "$meta.id"
    label 'process_single'

    conda "${projectDir}/modules/local/seqscreen/validate_fasta/environment.yml"

    input:
    tuple val(meta), path(fasta)
    path assets

    output:
    tuple val(meta), path("*.validated.fasta"), emit: fasta
    tuple val(meta), path("seqscreen.log")    , emit: log
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p logs
    echo -n " #### Launching SeqScreen DSL2 pipeline ....... " | tee seqscreen.log
    date '+%H:%M:%S %Y-%m-%d' | tee -a seqscreen.log

    python3 ${assets}/scripts/validate_fasta.py \\
        -f ${fasta} \\
        --max_seq_size=1000000000

    ln -s ${fasta} ${prefix}.validated.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqscreen: \$(seqscreen --version 2>&1 | sed 's/^SeqScreen v//; s/^.* //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.validated.fasta
    touch seqscreen.log
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqscreen: "stub"
    END_VERSIONS
    """
}
