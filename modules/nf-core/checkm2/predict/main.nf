// custom somatem module
process CHECKM2_PREDICT {
    tag "${meta.id}"
    label 'process_medium'

    // Note: Publishing is now handled globally via output blocks in main.nf
    // This provides centralized control and better organization of outputs

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(fasta, stageAs: "input_bins/*")
    path(db)

    output:
    tuple val(meta), path("${prefix}")                   , emit: checkm2_output
    tuple val(meta), path("${prefix}_checkm2_report.tsv"), emit: checkm2_tsv
    path("versions.yml")                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    checkm2 \\
        predict \\
        --input ${fasta} \\
        --output-directory ${prefix} \\
        --threads ${task.cpus} \\
        --database_path ${db} \\
        ${args}

    cp ${prefix}/quality_report.tsv ${prefix}_checkm2_report.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        checkm2: \$(checkm2 --version)
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}/
    touch ${prefix}_checkm2_report.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        checkm2: \$(checkm2 --version)
    END_VERSIONS
    """
}