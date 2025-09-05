// nf-core module file template copied from gms-16S/emu module
// refer there for more comments and stuff I deleted to maintain brevity..


process LEMUR {
    label 'process_high'
    
    conda "bioconda::lemur" // peg version with bioconda::lemur=1.0.1

    // optional: More reproducible than conda
    container "oras://community.wave.seqera.io/library/lemur:1.0.1--8e0c5d342d286d0b" 

    input:
      tuple val(meta), path(reads)

    output:
      tuple val(meta), path(output_dir)                           , emit: dir
      tuple val(meta), path("*relative_abundance.tsv")            , emit: report
      path "versions.yml"                                         , emit: versions

    script:
    def args = task.ext.args ?: ''
    output_dir = "${reads.baseName}-lemur-output"
    
    """
    lemur -i ${reads} \
      ${args} \
      --num-threads $task.cpus \
      -o ${output_dir}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lemur: \$(echo \$(lemur --version 2>&1) | sed 's/^.*lemur //; s/Using.*\$//' )
    END_VERSIONS
    """
}

