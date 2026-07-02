process SEQSCREEN_DIAMOND_BLASTX {
    tag "$meta.id"
    label 'process_high'

    conda "${projectDir}/modules/local/seqscreen/diamond_blastx/environment.yml"

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
    def evalue = params.seqscreen_evalue ?: 30
    """
    diamond blastx \\
        -q ${fasta} \\
        -d ${db}/diamond/uniref.mini.dmnd \\
        -o ${prefix}.daa \\
        --evalue ${evalue} \\
        --threads ${task.cpus} \\
        --block-size 200 \\
        --index-chunks 1 \\
        --salltitles \\
        --more-sensitive \\
        --min-orf 10 \\
        --masking 0 \\
        --top 5 \\
        -f 100

    diamond view \\
        -a ${prefix}.daa \\
        --top 5 \\
        --out ${prefix}.ur100.xml \\
        --outfmt 5

    diamond view \\
        -a ${prefix}.daa \\
        --top 5 \\
        --out ${prefix}.ur100.btab \\
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore ppos qframe score salltitles

    ln -s ${prefix}.ur100.btab functional_link.ur100.btab
    ln -s ${prefix}.ur100.xml functional_link.ur100.xml

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        diamond: \$(diamond version 2>&1 | head -n 1 | sed 's/^diamond version //')
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
        diamond: "stub"
    END_VERSIONS
    """
}
