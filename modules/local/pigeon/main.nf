process PIGEON {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(gfa), path(assembly), path(bins, stageAs: "bins/*")

    output:
    tuple val(meta), path("${prefix}/report.html")        , emit: html
    tuple val(meta), path("${prefix}/novel_metrics.json") , emit: metrics
    tuple val(meta), path("${prefix}")                    , emit: outdir
    path "versions.yml"                                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}_pigeon"
    
    """
    # List the bins directory to verify files are there
    echo "Contents of bins directory:"
    ls -la bins/
    
    # Check if gfatools is available
    if ! command -v gfatools &> /dev/null; then
        echo "ERROR: gfatools is not available. Please install gfatools."
        exit 1
    fi
    
    # Run pigeon analysis using GFA file instead of reads
    python ${projectDir}/bin/pigeon.py \\
        --gfa ${gfa} \\
        --assembly ${assembly} \\
        --bins_dir ${bins} \\
        --outdir ${prefix} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pigeon: "1.0.0"
        sourmash: \$(sourmash --version 2>&1 | sed 's/sourmash //')
        gfatools: \$(gfatools 2>&1 | head -n1 | sed 's/Version: //')
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}_pigeon"
    """
    echo $args  
    
    mkdir -p ${prefix}
    touch ${prefix}/report.html
    touch ${prefix}/novel_metrics.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pigeon: "1.0.0"
        sourmash: \$(sourmash --version 2>&1 | sed 's/sourmash //' || echo "N/A")
        gfatools: \$(gfatools 2>&1 | head -n1 | sed 's/Version: //' || echo "N/A")
        python: \$(python --version 2>&1 | sed 's/Python //' || echo "N/A")
    END_VERSIONS
    """
}