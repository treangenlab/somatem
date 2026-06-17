process SEQSCREEN_BOWTIE2 {
    tag "$meta.id:$target"
    label 'process_medium'

    conda "${projectDir}/modules/local/seqscreen/environment.yml"

    input:
    tuple val(meta), path(fasta)
    path db
    path assets
    val target

    output:
    tuple val(meta), path("*.sam"), emit: sam
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def db_path = target == 'vfdb' ? "${db}/bowtie2/vfdb/vfdb" : "${db}/bowtie2/bsat_ccl/blacklist.seqs.nt"
    def out = target == 'vfdb' ? "blacklist_vfdb.sam" : "blacklist_bsat.sam"
    """
    ${assets}/modules/bowtie2.sh \\
        --fasta=${fasta} \\
        --database=${db_path} \\
        --out=${out} \\
        --threads=${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(bowtie2 --version 2>&1 | head -n 1 | sed 's/^.*version //')
    END_VERSIONS
    """

    stub:
    """
    touch blacklist_${target}.sam
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: "stub"
    END_VERSIONS
    """
}
