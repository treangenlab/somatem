process SINGLEM_PIPE {
    tag "${meta.id}_${sample_type}"
    label 'process_medium'

    // Note: Publishing is now handled globally via output blocks in main.nf
    // This provides centralized control and better organization of outputs

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(reads)
    path metapackage
    val sample_type 

    output:
    tuple val(meta), path("*.tsv")                                      , emit: taxonomic_profile, optional: true
    tuple val(meta), path("*.html")                                     , emit: krona_profile, optional: true
    tuple val(meta), path("*_otu_table.csv")                           , emit: otu_table, optional: true
    tuple val(meta), path("*_archive_otu_table.csv")                   , emit: archive_otu_table, optional: true
    tuple val(meta), path("*.jplace")                                   , emit: jplace, optional: true
    path "versions.yml"                                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_${sample_type}"
    
    // Handle different input types
    def input_args = ""
    if (meta.single_end) {
        input_args = "--sequences ${reads}"
    } else {
        def (forward, reverse) = reads
        input_args = "--forward ${forward} --reverse ${reverse}"
    }
    
    // Handle genome input
    if (sample_type == 'genome') {
        input_args = "--genome-fasta-files ${reads.join(' ')}"
    }
    
    // Handle SRA input
    if (sample_type == 'sra') {
        input_args = "--sra-files ${reads.join(' ')}"
    }
    
    // Build metapackage argument (database)
    def metapackage_arg = metapackage ? "--metapackage ${metapackage}" : ""
    
    // Build output arguments
    def output_args = ""
    def has_output = false
    
    if (args.contains('--taxonomic-profile') || task.ext.taxonomic_profile) {
        output_args += " --taxonomic-profile ${prefix}_profile.tsv"
        has_output = true
    }
    if (args.contains('--taxonomic-profile-krona') || task.ext.krona_profile) {
        output_args += " --taxonomic-profile-krona ${prefix}_krona.html"
        has_output = true
    }
    if (args.contains('--otu-table') || task.ext.otu_table) {
        output_args += " --otu-table ${prefix}_otu_table.csv"
        has_output = true
    }
    if (args.contains('--archive-otu-table') || task.ext.archive_otu_table) {
        output_args += " --archive-otu-table ${prefix}_archive_otu_table.csv"
        has_output = true
    }
    if (args.contains('--output-jplace') || task.ext.jplace) {
        output_args += " --output-jplace ${prefix}"
        has_output = true
    }
    
    // If no outputs specified, provide defaults based on input type
    if (!has_output) {
        if (sample_type == 'reads') {
            output_args += " --taxonomic-profile ${prefix}_profile.tsv --otu-table ${prefix}_otu_table.csv"
        } else {
            output_args += " --otu-table ${prefix}_otu_table.csv"
        }
    }

    """
    singlem pipe \\
        ${input_args} \\
        ${metapackage_arg} \\
        ${output_args} \\
        --threads ${task.cpus} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        singlem: \$(singlem --version 2>&1 | sed 's/singlem //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_${sample_type}"
    """
    # Create default outputs based on input type if none specified
    has_output=false
    
    if [[ "${args}" == *"--taxonomic-profile"* ]] || [[ "${task.ext.taxonomic_profile}" == "true" ]]; then
        touch ${prefix}_profile.tsv
        has_output=true
    fi
    if [[ "${args}" == *"--taxonomic-profile-krona"* ]] || [[ "${task.ext.krona_profile}" == "true" ]]; then
        touch ${prefix}_krona.html
        has_output=true
    fi
    if [[ "${args}" == *"--otu-table"* ]] || [[ "${task.ext.otu_table}" == "true" ]]; then
        touch ${prefix}_otu_table.csv
        has_output=true
    fi
    if [[ "${args}" == *"--archive-otu-table"* ]] || [[ "${task.ext.archive_otu_table}" == "true" ]]; then
        touch ${prefix}_archive_otu_table.csv
        has_output=true
    fi
    if [[ "${args}" == *"--output-jplace"* ]] || [[ "${task.ext.jplace}" == "true" ]]; then
        touch ${prefix}.jplace
        has_output=true
    fi
    
    # If no outputs specified, create defaults
    if [[ "\$has_output" == "false" ]]; then
        if [[ "${sample_type}" == "reads" ]]; then
            touch ${prefix}_profile.tsv
            touch ${prefix}_otu_table.csv
        else
            touch ${prefix}_otu_table.csv
        fi
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        singlem: \$(singlem --version 2>&1 | sed 's/singlem //' || echo "0.19.0")
    END_VERSIONS
    """
}