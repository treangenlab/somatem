process TAXBURST {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/taxburst:0.3.0--pyhdfd78af_0':
        'biocontainers/taxburst:0.3.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(classification_file)
    val input_format

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.json"), emit: json, optional: true
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def format_arg = input_format ? "--input-format ${input_format}" : ""
    def save_json = args.contains('--save-json') ? "--save-json ${prefix}_taxburst.json" : ""
    
    """
    taxburst \\
        ${format_arg} \\
        ${args} \\
        ${save_json} \\
        -o ${prefix}_taxburst.html \\
        ${classification_file}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        taxburst: \$(micromamba list taxburst | grep -E 'taxburst\s+' | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_taxburst.html
    touch ${prefix}_taxburst.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        taxburst: \$(micromamba list taxburst | grep -E 'taxburst\s+' | awk '{print \$2}')
    END_VERSIONS
    """
}