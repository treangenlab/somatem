process VIRSORTER2 {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::virsorter=2.2.4"

    input:
    tuple val(meta), path(contigs)
    val min_len

    output:
    tuple val(meta), path("vs2_output/final-viral-combined.fa"), emit: viral
    tuple val(meta), path("vs2_output/final-viral-score.tsv")  , emit: scores
    tuple val(meta), path("vs2_output/final-viral-boundary.tsv"), emit: boundaries
    path "versions.yml"                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    virsorter run \\
        -w vs2_output \\
        -i ${contigs} \\
        --min-length ${min_len} \\
        -j ${task.cpus} \\
        --include-groups dsDNAphage,NCLDV,RNA,ssDNA,lavidaviridae \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        virsorter2: \$(virsorter --version 2>&1 | sed 's/VirSorter //g')
    END_VERSIONS
    """
}