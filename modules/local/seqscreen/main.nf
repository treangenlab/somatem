process SEQSCREEN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'quay.io/biocontainers/YOUR-TOOL-HERE' }"

    input:
    tuple val(meta), path(fasta)
    path(db)

    output:
    tuple val(meta), path("seqscreen_output"), emit: output_dir
    tuple val(meta), path("seqscreen_output/report_generation/seqscreen_report.tsv"), emit: report
    tuple val(meta), path("seqscreen_output/taxonomic_identification/taxonomic_assignment/taxonomic_results.txt"), emit: taxonomic_results
    tuple val(meta), path("seqscreen_output/seqscreen.log"), emit: log
    
    tuple val("${task.process}"), val('seqscreen'), eval("seqscreen --version"), topic: versions, emit: versions_seqscreen

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mode = params.seqscreen_mode ? " --${params.seqscreen_mode}" : ' --fast'

    """
    seqscreen \\
        $args \\
        --threads $task.cpus \\
        --fasta ${fasta} \\
        ${mode} \\
        --databases ${db} \\
        --working seqscreen_output
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo $args
    
    mkdir -p seqscreen_output/functional_annotation/functional_assignments
    mkdir -p seqscreen_output/taxonomic_identification/taxonomic_assignment
    mkdir -p seqscreen_output/report_generation
    
    touch seqscreen_output/report_generation/seqscreen_report.tsv
    touch seqscreen_output/functional_annotation/functional_assignments/functional_results.txt
    touch seqscreen_output/taxonomic_identification/taxonomic_assignment/taxonomic_results.txt
    touch seqscreen_output/seqscreen.log
    """
}
