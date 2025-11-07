process PIGEON {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(gfa), path(assembly), path(bins, stageAs: "bins/")

    output:
    tuple val(meta), path("${prefix}_report.html")     , emit: report
    tuple val(meta), path("${prefix}_metrics.json")    , emit: metrics
    path "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    echo "Contents of bins directory:"
    ls -la ${bins}

    # Run pigeon analysis
    pigeon.py \\
        --gfa ${gfa} \\
        --assembly ${assembly} \\
        --bins_dir ${bins} \\
        --outdir ./ \\
        $args

    # Append prefix to outputs
    cp "report.html" "${prefix}_report.html"
    cp "novel_metrics.json" "${prefix}_metrics.json"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pigeon: "1.0.0"
        sourmash: \$(sourmash --version 2>&1 | sed 's/sourmash //' | head -1)
        gfatools: \$(gfatools 2>&1 | head -n1 | sed 's/Version: //' || echo "unknown")
        python: \$(python --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo $args

    mkdir -p ${prefix}
    touch ${prefix}_report.html
    touch ${prefix}_metrics.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pigeon: "1.0.0"
        sourmash: "4.8.6"
        gfatools: "0.5"
        python: "3.11.0"
    END_VERSIONS
    """
}