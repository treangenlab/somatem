process SEQSCREEN_REPORT_HTML {
    tag "$meta.id"
    label 'process_medium'

    conda "${projectDir}/modules/local/seqscreen/environment.yml"

    input:
    tuple val(meta), path(fasta), path(report), path(blastx_xml)
    path db
    path assets
    val mode
    val rflag

    output:
    tuple val(meta), path("seqscreen_html_report"), emit: html
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    ${assets}/modules/html_report_generation.sh \\
        --report=${report} \\
        --fasta=${fasta} \\
        --blastx=${blastx_xml} \\
        --version=4.4 \\
        --mode=${mode} \\
        --rflag=${rflag} \\
        --gonetwork=${db}/go/go_network.txt \\
        --out=seqscreen_html_report/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqscreen: \$(seqscreen --version 2>&1 | sed 's/^SeqScreen v//; s/^.* //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p seqscreen_html_report
    touch seqscreen_html_report/index.html
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqscreen: "stub"
    END_VERSIONS
    """
}
