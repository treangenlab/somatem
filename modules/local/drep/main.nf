// dRep module for genome dereplication
// De-replicates a set of genomes by identifying highly similar genomes
// and selecting the best representative from each cluster

process DREP_DEREPLICATE {
    tag "$meta.id"
    label 'process_high'
    
    conda "bioconda::drep"

    input:
    tuple val(meta), path(genomes) // Multiple genome files to dereplicate
    val(ani_threshold)              // ANI threshold for clustering (default: 0.95)

    output:
    tuple val(meta), path("drep_output/dereplicated_genomes/*.fa*"), emit: genomes
    tuple val(meta), path("drep_output/data_tables/*.csv")         , emit: tables
    tuple val(meta), path("drep_output/figures/*.pdf")             , emit: figures, optional: true
    tuple val(meta), path("drep_output/log/*.log")                 , emit: logs
    path "versions.yml"                                            , emit: versions

    script:
    def args = task.ext.args ?: ''
    def ani = ani_threshold ?: 0.95
    
    """
    # Create input directory for genomes
    mkdir -p input_genomes
    
    # Copy/link all input genomes to input directory
    for genome in ${genomes}; do
        cp "\$genome" input_genomes/
    done
    
    # Run dRep dereplicate
    dRep dereplicate drep_output \\
        -g input_genomes/*.fa* \\
        -p ${task.cpus} \\
        -sa ${ani} \\
        --ignoreGenomeQuality \\
        ${args}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        drep: \$(dRep --version 2>&1 | sed 's/dRep //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p drep_output/dereplicated_genomes
    mkdir -p drep_output/data_tables
    mkdir -p drep_output/log
    touch drep_output/dereplicated_genomes/${prefix}_representative.fasta
    touch drep_output/data_tables/Cdb.csv
    touch drep_output/data_tables/Wdb.csv
    touch drep_output/log/logger.log
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        drep: \$(dRep --version 2>&1 | sed 's/dRep //')
    END_VERSIONS
    """
}
