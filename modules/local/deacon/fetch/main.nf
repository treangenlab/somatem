process DEACON_FETCH {
    tag "${index_name}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "quay.io/biocontainers/deacon:0.15.0--hdd79491_0"

    input:
    val index_name

    output:
    path "${index_name}.k31w15.idx", emit: index
    path 'versions.yml', emit: versions

    script:
    """
    deacon index fetch ${index_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deacon: \$(deacon --version | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    """
    touch ${index_name}.k31w15.idx
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deacon: "stub"
    END_VERSIONS
    """
}
