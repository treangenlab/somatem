/*
 * CLAIR3 MODULE
 * Calls variants from aligned reads using Clair3
 * Based on marlin.py call_variants_clair3() function
 */

process CLAIR3 {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(bam), path(bai), path(reference)
    val clair3_model
    val clair3_platform

    output:
    tuple val(meta), path("*.vcf.gz"), path("*.vcf.gz.tbi"), emit: vcf
    path "versions.yml", emit: versions

    script:
    def prefix = "${meta.id}"
    def model_arg = clair3_model ? "--model_path=${clair3_model}" : ""
    def model_path = clair3_model && clair3_model.name != 'input' ? clair3_model : "\${CONDA_PREFIX}/bin/models/${clair3_platform}"
    """
    # Ensure reference FASTA index exists
    if [ ! -f "${reference}.fai" ]; then
        samtools faidx ${reference}
    fi

    # Run Clair3 variant calling
    run_clair3.sh \\
        --bam_fn=${bam} \\
        --ref_fn=${reference} \\
        --output=clair3_output \\
        --threads=${task.cpus} \\
        --platform=${clair3_platform} \\
        ${model_arg}

    # Convert to uncompressed VCF, then re-compress and index
    bcftools view -Ov clair3_output/merge_output.vcf.gz > ${prefix}.vcf
    
    # Sort, bgzip, and tabix index
    bcftools sort -Oz -o ${prefix}.vcf.gz ${prefix}.vcf
    tabix -p vcf ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clair3: \$(run_clair3.sh --version 2>&1 | head -n1 || echo "unknown")
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = "${meta.id}"
    """
    touch ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clair3: stub
    END_VERSIONS
    """
}
