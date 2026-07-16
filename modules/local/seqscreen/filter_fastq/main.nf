process SEQSCREEN_FILTER_FASTQ {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/chopper:0.13.0--h7f49ad2_0' :
        'quay.io/biocontainers/chopper:0.13.0--h7f49ad2_0' }"

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("*.q15.fastq"), emit: fastq
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    chopper \\
        --input ${fastq} \\
        --threads $task.cpus \\
        --quality 15 \\
        > ${prefix}.q15.fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        chopper: \$(chopper --version 2>&1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.q15.fastq
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        chopper: "stub"
    END_VERSIONS
    """
}
