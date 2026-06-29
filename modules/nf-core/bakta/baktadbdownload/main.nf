// updated from nf-core module: update bakta to version 1.12.0
process BAKTA_BAKTADBDOWNLOAD {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bakta:1.12.0--pyhdfd78af_0' :
        'biocontainers/bakta:1.12.0--pyhdfd78af_0' }"
        // container for bakta v1.12.0 not tested! 

    output:
    path "db*"              , emit: db
    path "versions.yml"     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def db_type = (args =~ /(^|\s)--type(\s|=|$)/).find() ? '' : '--type light'
    """
    bakta_db \\
        download \\
        $db_type \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bakta: \$(echo \$(bakta_db --version) 2>&1 | cut -f '2' -d ' ')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def db_type = (args =~ /(^|\s)--type(\s|=|$)/).find() ? '' : '--type light'
    """
    echo "bakta_db \\
        download \\
        $db_type \\
        $args"

    mkdir db

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bakta: \$(echo \$(bakta_db --version) 2>&1 | cut -f '2' -d ' ')
    END_VERSIONS
    """
}
