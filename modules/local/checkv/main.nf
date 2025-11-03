process CHECKV {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::checkv=1.0.3"

    input:
    tuple val(meta), path(viral_contigs)

    output:
    tuple val(meta), path("checkv_output/proviruses.fna")     , emit: proviruses
    tuple val(meta), path("checkv_output/viruses.fna")       , emit: viruses
    tuple val(meta), path("checkv_output/quality_summary.tsv"), emit: quality
    tuple val(meta), path("checkv_output/curated_contigs.fna"), emit: curated
    path "versions.yml"                                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    checkv end_to_end \\
        ${viral_contigs} \\
        checkv_output \\
        -t ${task.cpus} \\
        ${args}

    # Combine viruses and proviruses into curated contigs
    cat checkv_output/viruses.fna checkv_output/proviruses.fna > checkv_output/curated_contigs.fna

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        checkv: \$(checkv -h 2>&1 | grep 'CheckV v' | sed 's/CheckV v//g')
    END_VERSIONS
    """
}