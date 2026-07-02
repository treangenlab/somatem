process SEQSCREEN_OUTLIER_DETECTION {
    tag "$meta.id"
    label 'process_single'

    conda "${projectDir}/modules/local/seqscreen/outlier_detection/environment.yml"

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
    def clean_btab = "${blastn_btab}_outlier_clean.btab"
    """
    sort \\
        -r \\
        -k 1,1 \\
        -k 12,12 \\
        ${blastn_btab} \\
        > ${blastn_btab}.srt

    python3 ${assets}/scripts/outlier_detection/score_blast.py \\
        -q ${fasta} \\
        -b ${blastn_btab}.srt \\
        -o outlier_detection.txt \\
        > outlier_detection.log 2>&1

    python3 ${assets}/scripts/outlier_detection/remove_outliers.py \\
        --outliers outlier_detection.txt \\
        --btab ${blastn_btab} \\
        --out ${clean_btab} \\
        >> outlier_detection.log 2>&1

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
