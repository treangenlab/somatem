

process SYLPH {
    conda "bioconda::sylph" // peg version with bioconda::name=version

    // optional: More reproducible than conda
    // container "oras://community.wave.seqera.io/library/name:version--build"  // generate with `wave containerize`

    input:
      path reads
      path reference 
      val threads

    output:
      path output

    script:
      output = "${reads.baseName}-sylph-profile.tsv"
    """
    sylph profile \
      ${reference} \
      ${reads} \
      -t ${threads} \
      -o ${output}
    """
}

