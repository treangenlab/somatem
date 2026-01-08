process SELECT_VARIANT_CALLER {
    tag "$meta.id"
    label 'process_single'
    
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(depth_file)
    val min_cov_for_clair3

    output:
    tuple val(meta), path("caller_choice.txt"), emit: choice
    tuple val(meta), path("coverage_stats.txt"), emit: stats
    path "versions.yml"                        , emit: versions

    script:
    """
    python $PWD/bin/select_variant_caller.py \\
        --depth ${depth_file} \\
        --min-cov ${min_cov_for_clair3} \\
        --output caller_choice.txt \\
        --stats coverage_stats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
