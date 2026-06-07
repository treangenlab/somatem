process BTYPER3 {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"

    publishDir { "${params.outdir}/isolate_analysis/${meta.id}/typing/btyper3" },
        mode: params.publish_dir_mode,
        pattern: "btyper3_out/**"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("btyper3_out"), emit: outdir
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = params.btyper3_args ?: ''
    """
    btyper3 \\
        -i ${fasta} \\
        -o btyper3_out \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        btyper3: \$(btyper3 --version 2>&1 | sed 's/^.*BTyper3 //; s/^btyper3 //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p btyper3_out
    touch btyper3_out/btyper3_stub_results.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        btyper3: "stub"
    END_VERSIONS
    """
}
