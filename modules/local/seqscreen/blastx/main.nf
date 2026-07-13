process SEQSCREEN_BLASTX {
    tag "$meta.id"
    label 'process_high'

    conda "${projectDir}/modules/local/seqscreen/blastx/environment.yml"

    input:
    tuple val(meta), path(fasta)
    path db

    output:
    tuple val(meta), path("${meta.id}.ur100.btab") , emit: btab
    tuple val(meta), path("${meta.id}.ur100.xml")  , emit: xml
    tuple val(meta), path("functional_link.ur100.btab"), emit: functional_btab
    tuple val(meta), path("functional_link.ur100.xml") , emit: functional_xml
    path "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = "${meta.id}"
    def evalue = params.seqscreen_evalue ?: 10
    """
    blastx \\
        -query ${fasta} \\
        -db ${db}/blast/UNIREF100.mini \\
        -out ${prefix}.ur100.asn \\
        -evalue ${evalue} \\
        -num_threads ${task.cpus} \\
        -max_target_seqs 500 \\
        -task blastx \\
        -seg no \\
        -outfmt 11 \\
        > blastx.log 2>&1

    blast_formatter \\
        -archive ${prefix}.ur100.asn \\
        -out ${prefix}.ur100.xml \\
        -outfmt 5 \\
        >> blastx.log 2>&1

    blast_formatter \\
        -archive ${prefix}.ur100.asn \\
        -out ${prefix}.ur100.btab \\
        -outfmt "6 std ppos qframe score salltitles" \\
        >> blastx.log 2>&1

    rm ${prefix}.ur100.asn

    ln -s ${prefix}.ur100.btab functional_link.ur100.btab
    ln -s ${prefix}.ur100.xml functional_link.ur100.xml

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: "\$(blastx -version 2>&1 | head -n 1 | sed 's/^blastx: //')"
    END_VERSIONS
    """

    stub:
    def prefix = "${meta.id}"
    """
    touch ${prefix}.ur100.btab
    touch ${prefix}.ur100.xml
    ln -s ${prefix}.ur100.btab functional_link.ur100.btab
    ln -s ${prefix}.ur100.xml functional_link.ur100.xml
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: "stub"
    END_VERSIONS
    """
}
