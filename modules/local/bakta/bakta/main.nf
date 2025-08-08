process BAKTA_BAKTA {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bakta:1.10.4--pyhdfd78af_0' :
        'biocontainers/bakta:1.10.4--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)
    path db
    path proteins
    path prodigal_tf

    output:
    tuple val(meta), path("${prefix}.embl")             , emit: embl
    tuple val(meta), path("${prefix}.faa")              , emit: faa
    tuple val(meta), path("${prefix}.ffn")              , emit: ffn
    tuple val(meta), path("${prefix}.fna")              , emit: fna
    tuple val(meta), path("${prefix}.gbff")             , emit: gbff
    tuple val(meta), path("${prefix}.gff3")             , emit: gff
    tuple val(meta), path("${prefix}.hypotheticals.tsv"), emit: hypotheticals_tsv
    tuple val(meta), path("${prefix}.hypotheticals.faa"), emit: hypotheticals_faa
    tuple val(meta), path("${prefix}.tsv")              , emit: tsv
    tuple val(meta), path("${prefix}.txt")              , emit: txt
    tuple val(meta), path("${prefix}.inference.tsv")    , emit: inference_tsv
    tuple val(meta), path("${prefix}.png"), optional: true, emit: png
    tuple val(meta), path("${prefix}.svg"), optional: true, emit: svg
    tuple val(meta), path("${prefix}.json")             , emit: json
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    def proteins_opt = proteins ? "--proteins ${proteins[0]}" : ""
    def prodigal_tf = prodigal_tf ? "--prodigal-tf ${prodigal_tf[0]}" : ""
    
    // Check if this bin should have plots based on completeness threshold from params
    def skip_plot = ""
    if (meta.completeness != null && meta.completeness < params.completeness_threshold) {
        skip_plot = "--skip-plot"
        println "Skipping plots for ${meta.id} (completeness: ${meta.completeness}% < ${params.completeness_threshold}%)"
    } else if (meta.completeness != null && meta.completeness >= params.completeness_threshold) {
        println "Generating plots for ${meta.id} (completeness: ${meta.completeness}% >= ${params.completeness_threshold}%)"
    } else {
        // If completeness is not available, skip plots by default
        skip_plot = "--skip-plot"
        println "Skipping plots for ${meta.id} (completeness unknown)"
    }
    
    """
    bakta \\
        $fasta \\
        $args \\
        --threads $task.cpus \\
        --prefix $prefix \\
        --meta \\
        --compliant \\
        $proteins_opt \\
        $prodigal_tf \\
        $skip_plot \\
        --db $db

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bakta: \$(echo \$(bakta --version) 2>&1 | cut -f '2' -d ' ')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.embl
    touch ${prefix}.faa
    touch ${prefix}.ffn
    touch ${prefix}.fna
    touch ${prefix}.gbff
    touch ${prefix}.gff3
    touch ${prefix}.hypotheticals.tsv
    touch ${prefix}.hypotheticals.faa
    touch ${prefix}.tsv
    touch ${prefix}.txt
    touch ${prefix}.inference.tsv
    touch ${prefix}.json
    
    # Only create PNG/SVG if completeness >= threshold
    if [[ "${meta.completeness}" != "null" ]] && (( \$(echo "${meta.completeness} >= ${params.completeness_threshold}" | bc -l) )); then
        touch ${prefix}.png
        touch ${prefix}.svg
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bakta: \$(echo \$(bakta --version) 2>&1 | cut -f '2' -d ' ')
    END_VERSIONS
    """
}