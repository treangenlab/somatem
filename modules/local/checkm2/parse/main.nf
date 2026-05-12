// custom somatem module
process CHECKM2_PARSE {
    tag "$meta.id"
    label 'process_low'

    // Note: Publishing is now handled globally via output blocks in main.nf
    // This provides centralized control and better organization of outputs

    conda "conda-forge::python=3.9 conda-forge::pandas=2.0.3"

    input:
    tuple val(meta), path(checkm2_tsv)
    tuple val(meta2), path(bins)

    output:
    tuple val(meta), path("bins_with_completeness.csv"), emit: completeness_map
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Run checkm2 parse script
    parse_checkm2.py ${checkm2_tsv}

    # Create versions file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
        pandas: \$(python3 -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """

    stub:
    """
    touch bins_with_completeness.csv
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
        pandas: \$(python3 -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """
}
