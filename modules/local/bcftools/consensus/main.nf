/*
 * BCFTOOLS_CONSENSUS MODULE
 * Generates consensus sequence from reference and high-confidence variants
 * Based on marlin.py consensus generation step
 */

process BCFTOOLS_CONSENSUS {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(reference), path(vcf), path(tbi)

    output:
    tuple val(meta), path("*_consensus.fasta"), emit: fasta
    path "versions.yml", emit: versions

    script:
    def prefix = "${meta.id}"
    """
    # Generate consensus sequence from reference and variants
    bcftools consensus -f ${reference} -o ${prefix}_consensus.fasta ${vcf}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //')
    END_VERSIONS
    """

    stub:
    def prefix = "${meta.id}"
    """
    touch ${prefix}_consensus.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: stub
    END_VERSIONS
    """
}
