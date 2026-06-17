process SEQSCREEN_OUTLIER_DETECTION {
    tag "$meta.id"
    label 'process_single'

    conda "${projectDir}/modules/local/seqscreen/environment.yml"

    input:
    tuple val(meta), path(fasta), path(blastn_btab)
    path assets

    output:
    tuple val(meta), path("*.nt.btab_outlier_clean.btab"), emit: clean_btab
    tuple val(meta), path("outlier_detection.txt")       , emit: outliers
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    ${assets}/modules/outlier_detection.sh \\
        --fasta=${fasta} \\
        --btab=${blastn_btab} \\
        --out=outlier_detection.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqscreen: \$(seqscreen --version 2>&1 | sed 's/^SeqScreen v//; s/^.* //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.nt.btab_outlier_clean.btab
    touch outlier_detection.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqscreen: "stub"
    END_VERSIONS
    """
}
