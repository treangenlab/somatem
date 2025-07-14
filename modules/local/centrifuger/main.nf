

process CENTRIFUGER_CLASSIFY {
    conda "bioconda::centrifuger" // peg version with bioconda::name=version

    // optional: More reproducible than conda
    // container "oras://community.wave.seqera.io/library/name:version--build"  // generate with `wave containerize`

    input:
      path reads
      path db_dir
      val threads

    output:
      path output_tsv

    script:
      output_tsv = "${reads.baseName}-centrifuger-output.tsv"
    """
    # find the index prefix (dig into the directory to find the index prefix of the 1st out of 4 cfr files)
    # source: nf-core/centrifuge: https://github.com/nf-core/modules/blob/master/modules/nf-core/centrifuge/centrifuge/main.nf#L43
    index_prefix=`find -L ${db_dir} -name "*.1.cfr" -not -name "._*"  | sed 's/\\.1.cfr\$//'`

    centrifuger -u ${reads} \
      -x \$index_prefix \
      -t ${threads} \
    > ${output_tsv}
    """
}

