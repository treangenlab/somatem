// custom somatem module
process SAMTOOLS_COVERAGE {
    tag "$meta.id"
    label 'process_single'

    // Note: Publishing is now handled globally via output blocks in main.nf
    // This provides centralized control and better organization of outputs

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(input), path(input_index)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fqi)

    output:
    tuple val(meta), path("*.coverage.txt"), emit: coverage
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference = fasta && !fasta.name.startsWith('OPTIONAL_FILE') ? "--reference ${fasta}" : ""
    def fqi_param = fqi && !fqi.name.startsWith('OPTIONAL_FILE') ? "--fqi ${fqi}" : ""
    """
    samtools \\
        coverage \\
        $args \\
        $reference \\
        $fqi_param \\
        -o ${prefix}.coverage.txt \\
        $input

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.coverage.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' )
    END_VERSIONS
    """
}
