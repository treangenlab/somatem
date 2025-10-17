process EMU_DOWNLOAD_DB {
    label 'process_single'

    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda "conda-forge::osfclient"
    
    output:
    tuple path ("species_taxid.fasta"), path ("taxonomy.tsv"), emit: emu_db_files
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    osf \\
        $args \\
        --project 56uf7 \\
        fetch osfstorage/emu-prebuilt/emu.tar

    tar -xvf emu.tar

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        osf: \$(echo \$(osf --version 2>&1) | grep 'osf version' | cut -f3 -d ' ')
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
    
    touch species_taxid.fasta
    touch taxonomy.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        osf: \$(echo \$(osf --version 2>&1) | grep 'osf version' | cut -f3 -d ' ')
    END_VERSIONS
    """
}


process EMU_STAGE_DB {
    label 'process_single'
    
    input:
    tuple path ("species_taxid.fasta"), path ("taxonomy.tsv")
    
    output:
    path "emu_db/", emit: emu_db
    
    script:
    """
    mkdir emu_db
    mv species_taxid.fasta emu_db/
    mv taxonomy.tsv emu_db/
    """
    
    stub:
    """
    mkdir emu_db
    mv species_taxid.fasta emu_db/
    mv taxonomy.tsv emu_db/
    """
}