/*
 * BCFTOOLS_CALL MODULE
 * Variant calling using bcftools mpileup + call pipeline
 * Suitable for low coverage (<8x) alignments
 */

process BCFTOOLS_CALL {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(bam), path(bai), path(reference)

    output:
    tuple val(meta), path("*.vcf.gz"), path("*.vcf.gz.tbi"), emit: vcf
    path "versions.yml", emit: versions

    script:
    def prefix = "${meta.id}"
    def avail_mem = task.memory ? "-m ${task.memory.toGiga()}G" : ""
    
    """
    # Run bcftools mpileup + call for low-coverage variant calling
    # This is preferred over Clair3 when coverage is <8x
    # -I flag skips indel calling (SNPs only for more reliable low-coverage calls)
    bcftools mpileup \\
        --threads ${task.cpus} \\
        --fasta-ref ${reference} \\
        -I \\
        ${bam} \\
        | bcftools call \\
            --threads ${task.cpus} \\
            --multiallelic-caller \\
            --variants-only \\
            --output-type z \\
            --output ${prefix}.vcf.gz
    
    # Index the output VCF
    bcftools index --tbi --threads ${task.cpus} ${prefix}.vcf.gz
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -n1 | sed 's/bcftools //g')
    END_VERSIONS
    """

    stub:
    def prefix = "${meta.id}"
    """
    touch ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: stub
    END_VERSIONS
    """
}
