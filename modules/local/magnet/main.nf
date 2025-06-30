

process MAGNET {
    conda "bioconda::magnet=version" // peg version with bioconda::name=version

    // optional: More reproducible than conda
    container "oras://community.wave.seqera.io/library/name:version--build"  // generate with `wave containerize`

    input:
      path reads
      path classification
      path output

    output:
      path output

    script:
    """
    magnet.py -i ${reads} \
      -o ${output} \
      -c ${classification}
    """
}


