process SEQSCREEN_BLASTN {
    tag "$meta.id:$target"
    label 'process_high'

    conda "${projectDir}/modules/local/seqscreen/blastn/environment.yml"

    input:
    tuple val(meta), path(fasta)
    path db
    val target

    output:
    tuple val(meta), path("*.btab"), emit: btab
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def db_path = target == 'megares' ? "${db}/megares/megares_full" : "${db}/blast/nt/nt"
    def out = target == 'megares' ? "${prefix}.megares.btab" : "${prefix}.nt.btab"
    def evalue = params.seqscreen_evalue ?: 10
    """
    blastn \\
        -query ${fasta} \\
        -db ${db_path} \\
        -out ${out} \\
        -evalue ${evalue} \\
        -num_threads ${task.cpus} \\
        -max_target_seqs 500 \\
        -dust no \\
        -outfmt "6 std qseq sseq staxid" \\
        > blastn.log 2>&1

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(blastn -version 2>&1 | head -n 1 | sed 's/^blastn: //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.${target}.btab
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: "stub"
    END_VERSIONS
    """
}
