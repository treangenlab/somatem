process SEQSCREEN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'quay.io/biocontainers/YOUR-TOOL-HERE' }"

    input:
    tuple val(meta), path(fasta)
    tuple val(meta), path(db)

    output:
    // TODO nf-core: Named file extensions MUST be emitted for ALL output channels
    tuple val(meta), path("seqscreen_dir"), emit: seqscreen_dir
    // TODO nf-core: List additional required output channels/values here
    // TODO nf-core: Update the command here to obtain the version number of the software used in this module
    // TODO nf-core: If multiple software packages are used in this module, all MUST be added here
    //               by copying the line below and replacing the current tool with the extra tool(s)
    tuple val("${task.process}"), val('seqscreen'), eval("seqscreen --version"), topic: versions, emit: versions_seqscreen

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mode = params.seqscreen_mode ? " --${params.seqscreen_mode}" : ' --ont'

    """
    seqscreen \\
        $args \\
        -@ $task.cpus \\
        --fasta ${fasta} \\
        ${mode} \\
        --databases ${db} \\
        -o ${prefix}.seqscreen_dir
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo $args
    
    touch ${prefix}.seqscreen_dir
    """
}
