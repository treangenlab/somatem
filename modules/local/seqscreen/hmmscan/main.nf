process SEQSCREEN_HMMSCAN {
    tag "$meta.id"
    label 'process_medium'

    conda "${projectDir}/modules/local/seqscreen/hmmscan/environment.yml"

    input:
    tuple val(meta), path(fasta)
    path db

    output:
    tuple val(meta), path("pfam_hmm"), emit: outdir
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p pfam_hmm
    esl-translate ${fasta} > pfam_hmm/${prefix}.translated.fasta

    hmmscan \\
        --cpu ${task.cpus} \\
        -o pfam_hmm/${prefix}.txt \\
        --tblout pfam_hmm/${prefix}.tab \\
        -E 1e-3 \\
        ${db}/hmmscan/Pfam-A.hmm \\
        pfam_hmm/${prefix}.translated.fasta \\
        > hmmscan.log 2>&1

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hmmer: "\$(hmmscan -h 2>&1 | grep -m1 '^# HMMER' | sed 's/^# HMMER //; s/ .*\$//')"
    END_VERSIONS
    """

    stub:
    """
    mkdir -p pfam_hmm
    touch pfam_hmm/${meta.id}.translated.fasta
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hmmer: "stub"
    END_VERSIONS
    """
}
