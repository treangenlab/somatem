process SUGGEST_MAGNET_CUTOFF {
    tag "$meta.id"
    label 'process_single'
    
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(lemur_table)
    val target_frac

    output:
    tuple val(meta), path("magnet_cutoff.txt"), emit: cutoff
    tuple val(meta), path("cutoff_stats.txt")  , emit: stats
    path "versions.yml"                        , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    python $PWD/bin/suggest_magnet_cutoff.py \\
        --input ${lemur_table} \\
        --target-frac ${target_frac} \\
        --output magnet_cutoff.txt \\
        --stats cutoff_stats.txt \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """
}
