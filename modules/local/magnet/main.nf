

process MAGNET {
    tag "$meta.id"
    label 'process_high'
    
    conda "${moduleDir}/dependencies.yml" // for locked env use: locked-spec-file.txt

    // optional: More reproducible than conda
    // container "oras://community.wave.seqera.io/library/name:version--build"  // generate with `wave containerize`

    input:
      tuple val(meta), path(reads)
      path(classification)

    output:
      tuple val(meta), path("*cluster_representative.csv")        , emit: report
      path "versions.yml"                                         , emit: versions

    script:
    def args = task.ext.args ?: ''
    
    """
    python ${moduleDir}/magnet-repo/magnet.py \
      ${args} \
      --threads $task.cpus \
      -i ${reads} \
      -c ${classification} \
      -o ./

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        magnet: default
    END_VERSIONS
    """

    
}


