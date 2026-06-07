process PYPOLCA {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"

    publishDir { "${params.outdir}/isolate_analysis/${meta.id}/polishing/pypolca" },
        mode: params.publish_dir_mode,
        pattern: "pypolca_out/**"
    publishDir { "${params.outdir}/isolate_analysis/${meta.id}/assembly" },
        mode: params.publish_dir_mode,
        pattern: "*.fasta"

    input:
    tuple val(meta), path(assembly), path(short_reads_1), path(short_reads_2)

    output:
    tuple val(meta), path("*.pypolca.fasta"), emit: assembly
    tuple val(meta), path("pypolca_out")    , emit: outdir
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = params.pypolca_args ?: '--careful'
    def prefix = meta.id
    """
    pypolca run \\
        -a ${assembly} \\
        -1 ${short_reads_1} \\
        -2 ${short_reads_2} \\
        -o pypolca_out \\
        -p pypolca \\
        -t ${task.cpus} \\
        ${args}

    cp pypolca_out/pypolca_corrected.fasta ${prefix}.pypolca.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pypolca: \$(pypolca -V 2>&1 | sed 's/^pypolca //')
    END_VERSIONS
    """
}
