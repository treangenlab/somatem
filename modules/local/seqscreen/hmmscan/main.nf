process SEQSCREEN_HMMSCAN {
    tag "$meta.id"
    label 'process_medium'

    conda "${projectDir}/modules/local/seqscreen/environment.yml"

    input:
    tuple val(meta), path(fasta)
    path db
    path assets

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
    ${assets}/modules/hmmscan.sh \\
        --fasta=pfam_hmm/${prefix}.translated.fasta \\
        --database=${db}/hmmscan/Pfam-A.hmm \\
        --out=pfam_hmm/${prefix} \\
        --threads=${task.cpus} \\
        --evalue=1e-3

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hmmer: \$(hmmscan -h 2>&1 | grep -m1 '^# HMMER' | sed 's/^# HMMER //; s/ .*\$//')
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
