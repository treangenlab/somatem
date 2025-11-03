process VRHYME {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::vrhyme=1.1.0"

    input:
    tuple val(meta), path(contigs)
    tuple val(meta2), path(coverage)

    output:
    tuple val(meta), path("vrhyme_output/vRhyme_best_bins.*.fasta"), emit: bins
    tuple val(meta), path("vrhyme_output/vRhyme_membership.tsv")   , emit: membership
    tuple val(meta), path("vrhyme_output/vRhyme_bin_info.tsv")     , emit: bin_info
    path "versions.yml"                                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    vRhyme \\
        -i ${contigs} \\
        -c ${coverage} \\
        -o vrhyme_output \\
        -t ${task.cpus} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vrhyme: \$(vRhyme --version 2>&1 | grep 'vRhyme v' | sed 's/vRhyme v//g')
    END_VERSIONS
    """
}