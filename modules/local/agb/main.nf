

process AGB {
    tag "$meta.id"
    label 'process_low'

    conda "almiheenko::agb" // peg version with bioconda::name=version
    

    // optional: More reproducible than conda
    // container "oras://community.wave.seqera.io/library/agb:version--build"  // generate with `wave containerize`

    input:
      tuple val(meta), path(assembly_dir)
 

    output:
      tuple val(meta), path("agb-output")

    script:
      // output_dir = "agb-output"
      assembler_name = "flye"
    """
    agb.py -i ${assembly_dir} -a ${assembler_name}
    """
}

