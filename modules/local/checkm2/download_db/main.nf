

process checkm2_DOWNLOAD_DB {
    tag 'checkm2_download_db'
    label 'process_single'

    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'biocontainers/YOUR-TOOL-HERE' }"

    input:
    path checkm2_db

    output:
    path checkm2_db, emit: checkm2_db
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    
    """
    if [[ ! -d $checkm2_db ]]; then
        mkdir -p $checkm2_db
    fi

    checkm2 database \\
        --download \\
        $args \\
        --path $checkm2_db

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        checkm2: \$(checkm2 --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    
    // TODO nf-core: A stub section should mimic the execution of the original module as best as possible
    //               Have a look at the following examples:
    //               Simple example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bcftools/annotate/main.nf#L47-L63
    //               Complex example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bedtools/split/main.nf#L38-L54
    // TODO nf-core: If the module doesn't use arguments ($args), you SHOULD remove:
    //               - The definition of args `def args = task.ext.args ?: ''` above.
    //               - The use of the variable in the script `echo $args ` below.
    """
    echo $args

    if [[ ! -d $checkm2_db ]]; then
        mkdir -p $checkm2_db
    fi
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tst: \$(tst --version)
    END_VERSIONS
    """
}
