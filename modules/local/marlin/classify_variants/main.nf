process CLASSIFY_VARIANTS {
    tag "$meta.id"
    label 'process_single'
    
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.17--haef29d1_0':
        'biocontainers/bcftools:1.17--haef29d1_0' }"

    input:
    tuple val(meta), path(vcf), path(tbi), path(depth)
    val af_high
    val af_med
    val high_cov_frac
    val gap_min_cov_frac
    val gap_min_cov_floor
    
    output:
    tuple val(meta), path("*_classified.vcf.gz"), path("*_classified.vcf.gz.tbi"), emit: classified_vcf
    tuple val(meta), path("*_highconf.vcf.gz"), path("*_highconf.vcf.gz.tbi"), emit: highconf_vcf
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def af_h = af_high ?: 0.9
    def af_m = af_med ?: 0.7
    def hc_frac = high_cov_frac ?: 0.5
    def gap_frac = gap_min_cov_frac ?: 0.1
    def gap_floor = gap_min_cov_floor ?: 2
    
    """
    python $PWD/bin/classify_variants.py \\
        --vcf ${vcf} \\
        --depth ${depth} \\
        --prefix ${prefix} \\
        --af-high ${af_h} \\
        --af-med ${af_m} \\
        --high-cov-frac ${hc_frac}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //')
    END_VERSIONS
    """
}
