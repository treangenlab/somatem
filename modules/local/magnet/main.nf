process MAGNET {
    tag "$meta.id"
    label 'process_high'
    
    // Use official MAGnet conda environment from GitHub (i just copied it from there)
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(reads)
    path(classification)

    output:
    tuple val(meta), path("cluster_representative.csv"), emit: report
    tuple val(meta), path("reference_genomes")         , emit: ref_dir
    tuple val(meta), path("*.bam")                      , emit: bam, optional: true
    path "versions.yml"                                 , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    // Infer abundance column index from classification file
    // MAGnet expects -a <column_index> for the abundance column
    // Lemur uses column 1 (0-indexed) for relative_abundance
    def abundance_col = task.ext.abundance_col ?: '1'
    
    // Default min-abundance if not specified
    def min_abundance = task.ext.min_abundance ?: '0.01'
    
    // MAGnet mode (ont, pb, or ilmn)
    def mode = task.ext.mode ?: 'ont'
    
    """
    # MAGnet expects to run in the output directory
    # Call magnet exactly as marlin.py does it
    magnet \\
        -c ${classification} \\
        -i ${reads} \\
        -o . \\
        -m ${mode} \\
        -a ${abundance_col} \\
        --min-abundance ${min_abundance} \\
        --threads ${task.cpus} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        magnet: \$(magnet --version 2>&1 | grep -oP 'v?\\d+\\.\\d+\\.\\d+' || echo "unknown")
    END_VERSIONS
    """
}