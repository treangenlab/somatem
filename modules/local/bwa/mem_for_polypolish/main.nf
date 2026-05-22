process BWA_MEM_FOR_POLYPOLISH {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"

    publishDir { "${params.outdir}/isolate_analysis/${meta.id}/polishing/polypolish" },
        mode: params.publish_dir_mode,
        pattern: "*.{sam,fasta}"

    input:
    tuple val(meta), path(assembly), path(short_reads_1), path(short_reads_2)

    output:
    tuple val(meta), path("*.for_polypolish.fasta"), path("*.alignments_1.sam"), path("*.alignments_2.sam"), emit: alignments
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = meta.id
    """
    cp ${assembly} ${prefix}.for_polypolish.fasta
    bwa index ${prefix}.for_polypolish.fasta

    bwa mem \\
        -t ${task.cpus} \\
        -a \\
        ${prefix}.for_polypolish.fasta \\
        ${short_reads_1} \\
        > ${prefix}.alignments_1.sam

    bwa mem \\
        -t ${task.cpus} \\
        -a \\
        ${prefix}.for_polypolish.fasta \\
        ${short_reads_2} \\
        > ${prefix}.alignments_2.sam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(bwa 2>&1 | sed -n 's/^Version: //p' | head -n 1)
    END_VERSIONS
    """
}
