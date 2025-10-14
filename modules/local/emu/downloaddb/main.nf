process EMU_DOWNLOADDB {
    label 'process_single'

    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda "conda-forge::osfclient"
    
    output:
    path "db", emit: emu_db
    tuple path ("db/species_taxid.fasta"), path ("db/taxonomy.tsv"), emit: species_and_taxonomy

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    osf \\
        $args \\
        --project 56uf7 \\
        fetch osfstorage/emu-prebuilt/emu.tar.gz

    tar -xzf emu.tar.gz
    mv emu db
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
    
    mkdir db
    touch db/species_taxid.fasta
    touch db/taxonomy.tsv

    """
}
