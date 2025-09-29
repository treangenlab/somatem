process RHEA {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'biocontainers/YOUR-TOOL-HERE' }"

    input:
    tuple val(meta), path(multi_reads) // list of multiple fastq.gz read files

    output:
    tuple val(meta), path("metaflye/")                    , emit: assembly_dir
    tuple val(meta), path("metaflye/assembly_graph.gfa")  , emit: assembly_graph
    tuple val(meta), path("Bandage_metadata.csv")         , emit: bandage_metadata
    tuple val(meta), path("bp_counts.tsv")                , emit: bp_counts
    tuple val(meta), path("edge_coverage-*.tsv")          , emit: edge_coverage
    tuple val(meta), path("node_coverage_norm.csv")       , emit: node_coverage_norm
    tuple val(meta), path("node_coverage.csv")            , emit: node_coverage
    tuple val(meta), path("rhea.log")                     , emit: rhea_log
    tuple val(meta), path("structural_variants-c0.tsv")   , emit: structural_variants
    tuple val(meta), path("*.gaf")                        , emit: gaf_files
    
    path("versions.yml")                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    python ${moduleDir}/rhea-repo/rhea.py \\
        $args \\
        --threads $task.cpus \\
        ${multi_reads.join(' ')}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        flye: \$(flye --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    echo $args
    
    touch Bandage_metadata.csv
    touch bp_counts.tsv
    touch edge_coverage-c0.tsv
    touch metaflye
    touch metaflye/assembly_graph.gfa
    touch node_coverage_norm.csv
    touch node_coverage.csv
    touch rhea.log
    touch structural_variants-c0.tsv
    touch sample1.gaf


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        flye: \$(flye --version)
    END_VERSIONS
    """
}
