process SEQSCREEN_REPORT_TSV {
    tag "$meta.id"
    label 'process_single'

    conda "${projectDir}/modules/local/seqscreen/report_tsv/environment.yml"

    input:
    tuple val(meta), path(fasta), path(taxonomy), path(functional), path(bsat), path(vfdb)
    path db
    path assets
    val mode

    output:
    tuple val(meta), path("*seqscreen_report.tsv"), emit: report
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = params.seqscreen_report_prefix ? "${fasta.baseName}_" : ''
    """
    python3 ${assets}/scripts/seqscreen_tsv_report.py \\
        --taxonomy=${taxonomy} \\
        --functional=${functional} \\
        --taxlookup=${db}/taxonomy/taxa_lookup.txt \\
        --funsocs=${db}/funsocs \\
        --fasta=${fasta} \\
        --mode=${mode} \\
        --bsat=${bsat} \\
        --vfdb=${vfdb} \\
        --out=${prefix}seqscreen_report.tsv \\
        --parenttax=${db}/tax_to_parent.pck \\
        --mergedtax=${db}/merged_taxa.pck

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqscreen: "4.5"
    END_VERSIONS
    """

    stub:
    """
    touch seqscreen_report.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqscreen: "stub"
    END_VERSIONS
    """
}
