process SEQSCREEN_RAPSEARCH2 {
    tag "$meta.id:$target"
    label 'process_medium'

    conda "${projectDir}/modules/local/seqscreen/environment.yml"

    input:
    tuple val(meta), path(fasta)
    path db
    path assets
    val target

    output:
    tuple val(meta), path("*.m8"), emit: m8
    path "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def db_path = target == 'vfdb' ? "${db}/rapsearch2/vfdb/vfdb.seqs.aa" : "${db}/rapsearch2/bsat_ccl/blacklist.seqs.aa"
    def out = target == 'vfdb' ? "blacklist_vfdb" : "blacklist_bsat"
    """
    ${assets}/modules/rapsearch2.sh \\
        --fasta=${fasta} \\
        --database=${db_path} \\
        --out=${out} \\
        --evalue=1e-9 \\
        --threads=${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rapsearch2: \$(rapsearch -h 2>&1 | head -n 1 | sed 's/^.*RAPSearch2 //; s/ .*\$//' || true)
    END_VERSIONS
    """

    stub:
    """
    touch blacklist_${target}.m8
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rapsearch2: "stub"
    END_VERSIONS
    """
}
