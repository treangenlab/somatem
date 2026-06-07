process POLYPOLISH_FILTER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"

    publishDir { "${params.outdir}/isolate_analysis/${meta.id}/polishing/polypolish" },
        mode: params.publish_dir_mode,
        pattern: "*.{sam,fasta}"

    input:
    tuple val(meta), path(assembly), path(alignments_1), path(alignments_2)

    output:
    tuple val(meta), path("*.for_polypolish.fasta"), path("*.filtered_1.sam"), path("*.filtered_2.sam"), emit: alignments
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def filter_args = params.polypolish_filter_args ?: ''
    def prefix = meta.id
    """
    polypolish filter \\
        --in1 ${alignments_1} \\
        --in2 ${alignments_2} \\
        --out1 ${prefix}.filtered_1.sam \\
        --out2 ${prefix}.filtered_2.sam \\
        ${filter_args}

    cp ${assembly} ${prefix}.for_polypolish.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        polypolish: \$(polypolish --version 2>&1 | sed 's/^polypolish //')
    END_VERSIONS
    """
}
