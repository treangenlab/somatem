process SYLPH_DOWNLOAD_DB {
    tag "sylph_${db_name}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    input:
    val db_url
    val db_name

    output:
    path "${db_name}", emit: db
    path 'versions.yml', emit: versions

    script:
    """
    aria2c \
        --continue=true \
        --out ${db_name} \
        ${db_url}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        aria2: \$(aria2c --version | awk 'NR == 1 { print \$3 }')
    END_VERSIONS
    """

    stub:
    """
    touch ${db_name}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        aria2: "stub"
    END_VERSIONS
    """
}
