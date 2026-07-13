

process MAGNET {
    tag "$meta.id"
    label 'process_high'

    publishDir { "${params.outdir}/taxonomic_profiling/${meta.id}/magnet" }, mode: params.publish_dir_mode, pattern: '*.csv'
    
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
    magnet \
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

