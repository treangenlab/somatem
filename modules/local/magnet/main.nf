

process MAGNET {
    conda "${moduleDir}/spec-file.txt" // peg version with bioconda::name=version

    // optional: More reproducible than conda
    // container "oras://community.wave.seqera.io/library/name:version--build"  // generate with `wave containerize`

    input:
      path reads
      path classification

    output:
      path output_dir

    script:
      output_dir = "${reads.baseName}-magnet-output"
    """
    python ${moduleDir}/magnet-repo/magnet.py \
      -i ${reads} \
      -c ${classification} \
      -o ${output_dir}
    """
}


