

process TOOL {
    conda "bioconda::name=version" // peg version with bioconda::name=version

    // optional: More reproducible than conda
    container "oras://community.wave.seqera.io/library/name:version--build"  // generate with `wave containerize`

    input:
      path reads
      path database_dir
      path taxonomy
      val rank

    output:
      path output_dir

    script:
      output_dir = "${reads.baseName}-name-output"
    """
    name -i ${reads} \
      -o ${output_dir} \
      -d ${database_dir} \
      --tax-path ${taxonomy} \
      -r ${rank}
    """
}

