process SINGLEM_APPRAISE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/singlem:0.19.0--pyhdfd78af_0':
        'biocontainers/singlem:0.19.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(metagenome_otu_tables), path(genome_otu_tables), path(assembly_otu_tables)
    path metapackage

    output:
    tuple val(meta), path("*.txt")                                      , emit: summary
    tuple val(meta), path("*_binned.csv")         , optional: true      , emit: binned_otu_table
    tuple val(meta), path("*_unbinned.csv")       , optional: true      , emit: unbinned_otu_table
    tuple val(meta), path("*_assembled.csv")      , optional: true      , emit: assembled_otu_table
    tuple val(meta), path("*_unaccounted.csv")    , optional: true      , emit: unaccounted_otu_table
    tuple val(meta), path("*.svg")                , optional: true      , emit: plot
    path "versions.yml"                                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    // Build input arguments
    def metagenome_args = metagenome_otu_tables ? "--metagenome-otu-tables ${metagenome_otu_tables.join(' ')}" : ""
    def genome_args = genome_otu_tables ? "--genome-otu-tables ${genome_otu_tables.join(' ')}" : ""
    def assembly_args = assembly_otu_tables ? "--assembly-otu-tables ${assembly_otu_tables.join(' ')}" : ""
    def metapackage_args = metapackage ? "--metapackage ${metapackage}" : ""
    
    // Build output arguments
    def output_args = ""
    if (args.contains('--output-binned-otu-table') || task.ext.output_binned) {
        output_args += " --output-binned-otu-table ${prefix}_binned.csv"
    }
    if (args.contains('--output-unbinned-otu-table') || task.ext.output_unbinned) {
        output_args += " --output-unbinned-otu-table ${prefix}_unbinned.csv"
    }
    if (args.contains('--output-assembled-otu-table') || task.ext.output_assembled) {
        output_args += " --output-assembled-otu-table ${prefix}_assembled.csv"
    }
    if (args.contains('--output-unaccounted-for-otu-table') || task.ext.output_unaccounted) {
        output_args += " --output-unaccounted-for-otu-table ${prefix}_unaccounted.csv"
    }
    if (args.contains('--plot') || task.ext.plot) {
        output_args += " --plot ${prefix}_plot.svg"
    }

    """
    singlem appraise \\
        ${metagenome_args} \\
        ${genome_args} \\
        ${assembly_args} \\
        ${metapackage_args} \\
        ${output_args} \\
        ${args} \\
        > ${prefix}_appraise_summary.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        singlem: \$(singlem --version 2>&1 | sed 's/singlem //')
    END_VERSIONS
    """
}