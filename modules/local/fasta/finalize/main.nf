process FASTA_FINALIZE {
    tag "$meta.id"
    label 'process_single'

    publishDir { "${params.outdir}/isolate_analysis/${meta.id}/assembly" },
        mode: params.publish_dir_mode,
        pattern: "*.fasta"

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("*.fasta"), emit: assembly
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = meta.id
    """
    cp ${assembly} ${prefix}.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cp: \$(cp --version 2>&1 | head -n 1 | sed 's/^cp (GNU coreutils) //')
    END_VERSIONS
    """
}
