process SYLPH_DOWNLOAD_DB {
    label 'process_single'
   
    input:
    val (target_taxonomy) // string: target taxonomy for the database
    val (db_picker) // string: path to the CSV file containing database versions

    output:
    path ("*.syldb"), emit: sylph_db
    // path "versions.yml", emit: versions // not emitted since I don't have write access to the existing db directory

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    # Execute the script and capture its output
    db_version=\$(get_sylph_db_version.sh "${target_taxonomy}" "${db_picker}" 2>/dev/null)

    # print db version
    echo "Database version: \${db_version}"
    
    wget \\
        http://faust.compbio.cs.cmu.edu/sylph-stuff/\${db_version}.syldb

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sylph_DB: \${db_version}
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
    
    touch *.syldb

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(echo \$(wget --version 2>&1) | grep 'wget version' | cut -f3 -d ' ')
    END_VERSIONS
    """
}
