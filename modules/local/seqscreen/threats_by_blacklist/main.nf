process SEQSCREEN_THREATS_BY_BLACKLIST {
    tag "$meta.id:$target"
    label 'process_single'

    conda "${projectDir}/modules/local/seqscreen/threats_by_blacklist/environment.yml"

    input:
    tuple val(meta), path(sam), path(m8)
    path assets
    val target

    output:
    tuple val(meta), path("threats_by_${target}.txt"), emit: threats
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    python3 ${assets}/scripts/threats_by_blacklist.py \\
        --rapsearch-results=${m8} \\
        --bowtie2-results=${sam} \\
        --out=threats_by_${target}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqscreen: "4.5"
    END_VERSIONS
    """

    stub:
    """
    touch threats_by_${target}.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqscreen: "stub"
    END_VERSIONS
    """
}
