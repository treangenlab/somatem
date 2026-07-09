process SEQSCREEN_REPORT_HTML {
    tag "$meta.id"
    label 'process_medium'

    conda "${projectDir}/modules/local/seqscreen/report_html/environment.yml"

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
    python3 ${assets}/scripts/html_report_generation/generateHtmlReport.py \\
        -f ${fasta} \\
        -r ${report} \\
        -b ${blastx_xml} \\
        -o seqscreen_html_report/ \\
        --rflag ${rflag} \\
        --version 4.4 \\
        -d ${assets}/scripts/html_report_generation/libs/ \\
        -t ${assets}/scripts/html_report_generation/data/template.html \\
        -g ${assets}/scripts/html_report_generation/data/go_names.txt \\
        -n ${db}/go/go_network.txt \\
        --funsocs ${assets}/scripts/html_report_generation/data/funsocs_description.txt \\
        --mode ${mode} \\
        --go_template ${assets}/scripts/html_report_generation/data/go_template.html \\
        > html_report_generation.log 2>&1 || {
            cat html_report_generation.log >&2
            exit 1
        }

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqscreen: "4.5"
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
