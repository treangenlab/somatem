process IPHOP {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::iphop=1.4.1"

    input:
    tuple val(meta), path(viral_contigs)

    output:
    tuple val(meta), path("iphop_output/Host_prediction_to_genus_m90.csv"), emit: genus_predictions
    tuple val(meta), path("iphop_output/Host_prediction_to_species_m90.csv"), emit: species_predictions
    tuple val(meta), path("iphop_output/Detailed_output_by_tool.tsv"), emit: detailed_output
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    iphop predict \\
        --fa_file ${viral_contigs} \\
        --out_dir iphop_output \\
        --num_threads ${task.cpus} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iphop: \$(iphop --version 2>&1 | grep 'iPHoP v' | sed 's/iPHoP v//g')
    END_VERSIONS
    """
}