process SEQSCREEN_TAXONOMIC_ASSIGNMENT_FAST {
    tag "$meta.id"
    label 'process_single'

    conda "${projectDir}/modules/local/seqscreen/taxonomic_assignment/environment.yml"

    input:
    tuple val(meta), path(diamond_btab), path(centrifuge_results)
    path assets

    output:
    tuple val(meta), path("taxonomic_results.txt"), emit: results
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def taxlimit = params.seqscreen_taxlimit ?: 25
    """
    python3 ${assets}/scripts/consolidatedtax_2_taxreport_fast.py \\
        --diamond=${diamond_btab} \\
        --centrifuge=${centrifuge_results} \\
        --out=taxonomic_results.txt \\
        --cutoff=1 \\
        --taxlimit=${taxlimit}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqscreen: \$(seqscreen --version 2>&1 | sed 's/^SeqScreen v//; s/^.* //')
    END_VERSIONS
    """

    stub:
    """
    touch taxonomic_results.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqscreen: "stub"
    END_VERSIONS
    """
}

process SEQSCREEN_TAXONOMIC_ASSIGNMENT_SENSITIVE {
    tag "$meta.id"
    label 'process_single'

    conda "${projectDir}/modules/local/seqscreen/taxonomic_assignment/environment.yml"

    input:
    tuple val(meta), path(blastx_btab), path(clean_blastn_btab)
    path assets

    output:
    tuple val(meta), path("taxonomic_results.txt"), emit: results
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    python3 ${assets}/scripts/btab_2_tax_report_sensitive.py \\
        --blastx=${blastx_btab} \\
        --blastn=${clean_blastn_btab} \\
        --out=taxonomic_results.txt \\
        --cutoff=1

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqscreen: \$(seqscreen --version 2>&1 | sed 's/^SeqScreen v//; s/^.* //')
    END_VERSIONS
    """

    stub:
    """
    touch taxonomic_results.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqscreen: "stub"
    END_VERSIONS
    """
}
