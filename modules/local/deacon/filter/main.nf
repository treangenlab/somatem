process DEACON_FILTER {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "quay.io/biocontainers/deacon:0.15.0--hdd79491_0"

    input:
    tuple val(meta), path(reads)
    path index

    output:
    tuple val(meta), path('*.deacon.clean*.fastq.gz'), emit: fastq
    tuple val(meta), path('*.deacon.summary.json'), emit: json
    path 'versions.yml', emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: meta.id
    def read_list = [reads].flatten()
    def outputs = meta.single_end
        ? "-o ${prefix}.deacon.clean.fastq.gz"
        : "-o ${prefix}.deacon.clean_1.fastq.gz -O ${prefix}.deacon.clean_2.fastq.gz"
    """
    deacon filter \
        --deplete \
        --threads ${task.cpus} \
        --summary ${prefix}.deacon.summary.json \
        ${args} \
        ${index} \
        ${read_list.join(' ')} \
        ${outputs}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deacon: \$(deacon --version | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id
    def output_commands = meta.single_end
        ? "touch ${prefix}.deacon.clean.fastq.gz"
        : "touch ${prefix}.deacon.clean_1.fastq.gz ${prefix}.deacon.clean_2.fastq.gz"
    """
    ${output_commands}
    echo '{}' > ${prefix}.deacon.summary.json
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deacon: "stub"
    END_VERSIONS
    """
}
