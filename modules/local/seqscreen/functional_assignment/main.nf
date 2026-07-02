process SEQSCREEN_FUNCTIONAL_ASSIGNMENT {
    tag "$meta.id"
    label 'process_medium'

    conda "${projectDir}/modules/local/seqscreen/functional_assignment/environment.yml"

    input:
    tuple val(meta), path(fasta), path(blastx_btab)
    path db
    path assets

    output:
    tuple val(meta), path("functional_results.txt")      , emit: results
    tuple val(meta), path("functional_assignments.log")  , emit: log
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def ancestral = params.seqscreen_ancestral ? '--ancestral' : ''
    def includecc = params.seqscreen_includecc ? '--include-cc' : ''
    def cutoff = params.seqscreen_bitscore ?: 5
    """
    python3 ${assets}/scripts/functional_report_generation.py \\
        --fasta=${fasta} \\
        --urbtab=${blastx_btab} \\
        --go=${db}/go/go_network.txt \\
        --out=functional_results.txt \\
        --annotation=${db}/annotation_scores.pck \\
        ${ancestral} \\
        ${includecc} \\
        --cutoff=${cutoff} > functional_assignments.log 2>&1

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqscreen: \$(seqscreen --version 2>&1 | sed 's/^SeqScreen v//; s/^.* //')
    END_VERSIONS
    """

    stub:
    """
    touch functional_results.txt
    touch functional_assignments.log
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqscreen: "stub"
    END_VERSIONS
    """
}
