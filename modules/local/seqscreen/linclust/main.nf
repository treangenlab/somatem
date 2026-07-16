process SEQSCREEN_LINCLUST {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mmseqs2:18.8cc5c--hd6d6fdc_0' :
        'quay.io/biocontainers/mmseqs2:18.8cc5c--hd6d6fdc_0' }"

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("*.representatives.fasta"), emit: fasta
    path "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p tmp

    mmseqs createdb ${fastq} reads_db
    mmseqs linclust \\
        reads_db \\
        clusters_db \\
        tmp \\
        --min-seq-id 0.99 \\
        --threads $task.cpus
    mmseqs createsubdb clusters_db reads_db representatives_db
    mmseqs convert2fasta representatives_db ${prefix}.representatives.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mmseqs2: \$(mmseqs version)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.representatives.fasta
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mmseqs2: "stub"
    END_VERSIONS
    """
}
