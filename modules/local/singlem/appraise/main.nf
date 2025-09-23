process SINGLEM_APPRAISE {
    tag "$meta.id"
    label 'process_medium'

    // Outputs
    publishDir "${params.output_dir}/appraise/${meta.id}", mode: 'copy', pattern: "*.tsv"
    publishDir "${params.output_dir}/appraise/${meta.id}", mode: 'copy', pattern: "*.svg"

    conda "${moduleDir}/environment.yml"

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
    
    // Build input arguments - handle both single files and lists
    def metagenome_args = ""
    if (metagenome_otu_tables) {
        if (metagenome_otu_tables instanceof List) {
            metagenome_args = "--metagenome-otu-tables ${metagenome_otu_tables.join(' ')}"
        } else {
            metagenome_args = "--metagenome-otu-tables ${metagenome_otu_tables}"
        }
    }
    
    def genome_args = ""
    if (genome_otu_tables) {
        if (genome_otu_tables instanceof List) {
            genome_args = "--genome-otu-tables ${genome_otu_tables.join(' ')}"
        } else {
            genome_args = "--genome-otu-tables ${genome_otu_tables}"
        }
    }
    
    def assembly_args = ""
    if (assembly_otu_tables && assembly_otu_tables.size() > 0) {
        if (assembly_otu_tables instanceof List) {
            assembly_args = "--assembly-otu-tables ${assembly_otu_tables.join(' ')}"
        } else {
            assembly_args = "--assembly-otu-tables ${assembly_otu_tables}"
        }
    }
    
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

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_appraise_summary.txt
    touch ${prefix}_binned.csv
    touch ${prefix}_unbinned.csv
    touch ${prefix}_assembled.csv
    touch ${prefix}_unaccounted.csv
    touch ${prefix}_plot.svg

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        singlem: \$(singlem --version 2>&1 | sed 's/singlem //')
    END_VERSIONS
    """
}