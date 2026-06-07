// Plot assembly graph. Currently set to flye assembler.

process AGB {
    tag "$meta.id"
    label 'process_low'

    publishDir { "${params.output_dir}/assembly_graph/${meta.id}" }, mode: 'copy', pattern: "agb_output/**"

    conda "almiheenko::agb" // peg version with bioconda::name=version
    

    // optional: More reproducible than conda
    // container "oras://community.wave.seqera.io/library/agb:version--build"  // generate with `wave containerize`

    input:
      tuple val(meta), path(gfa)
      tuple val(_meta1), path(gv)
      tuple val(_meta2), path(txt) 

    output:
      tuple val(meta), path("agb_output/"), emit: assembly_graph

    script:
      assembler_name = "flye"
    """
    # stage the gfa, gv and txt files into a directory
    mkdir -p agb_input
    gunzip -c ${gfa} > agb_input/assembly_graph.gfa
    gunzip -c ${gv} > agb_input/assembly_graph.gv
    mv ${txt} agb_input/assembly_info.txt
    
    agb.py \\
    -a ${assembler_name} \\
    -i agb_input \\
    -t $task.cpus
    """
}

