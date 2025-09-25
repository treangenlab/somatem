process TAXBURST_CONVERT {
    tag   "$meta.id"
    label 'process_single'

    conda      "${moduleDir}/environment.yml"
    container  "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
                   'https://depot.galaxyproject.org/singularity/python:3.9--h2a3dd35_0' :
                   'quay.io/biocontainers/python:3.9--h2a3dd35_0' }"

    input:
    tuple val(meta), path(classification_file)
    val tool_name
    // for centrifuger mode we need the taxonomy DB:
    // val tax_db_dir        // e.g. params.centrifuger_db

    output:
    tuple val(meta), path("${meta.id}.krona.tsv"), emit: converted
    path "versions.yml", emit: versions

    /*
     * ext.prefix can override the output basename,
     * else we use meta.id
     */
    script:
    def prefix = task.ext.prefix ?: meta.id

    """
    # invoke converter script
    taxburst_prep.py \\
      --from_tool ${tool_name} \\
      --format krona \\
      --output ${prefix}.krona.tsv \\
      ${classification_file}

    # record tool versions
    echo "${task.process}:" > versions.yml
    echo "  python: \$(python --version | sed 's/Python //'\)" >> versions.yml
    """
}
