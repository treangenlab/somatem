// nf-core module file template copied from gms-16S/emu module
// refer there for more comments and stuff I deleted to maintain brevity..


process lemur {
    conda "bioconda::lemur" // issue with bioconda::lemur:1.0.1

    // optional: More reproducible than conda
    container "oras://community.wave.seqera.io/library/lemur:1.0.1--8e0c5d342d286d0b" 
    
    input:
      path reads
      path database_dir
      path taxonomy
      val rank

    output:
      path output_dir

    script:
      output_dir = "${reads.baseName}-lemur-output"
    """
    lemur -i ${reads} \
      -o ${output_dir} \
      -d ${database_dir} \
      --tax-path ${taxonomy} \
      -r ${rank}
    """
}

