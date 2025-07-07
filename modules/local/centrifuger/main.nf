

process CENTRIFUGER_CLASSIFY {
    conda "bioconda::centrifuger" // peg version with bioconda::name=version

    // optional: More reproducible than conda
    // container "oras://community.wave.seqera.io/library/name:version--build"  // generate with `wave containerize`

    input:
      path reads
      path index_prefix
      val threads

    output:
      path output_tsv

    script:
      output_tsv = "${reads.baseName}-centrifuger-output.tsv"
    """
    centrifuger -u ${reads} \
      -x ${index_prefix} \
      -t ${threads} \
    > ${output_tsv}
    """
}

