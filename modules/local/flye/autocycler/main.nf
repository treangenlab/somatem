process FLYE_AUTOCYCLER {
    tag "${meta.id}:${meta.subsample_id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"

    publishDir { "${params.outdir}/isolate_analysis/${meta.id}/autocycler/assemblies" },
        mode: params.publish_dir_mode,
        pattern: "*.{fasta,gfa,log,txt}"

    input:
    tuple val(meta), path(subsample_reads), path(genome_size_file), val(species), val(read_type)

    output:
    tuple val(meta), path("*.assembly.fasta"), emit: assembly
    tuple val(meta), path("*.assembly.gfa")  , optional: true, emit: graph
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def mode = params.autocycler_flye_mode ?: (!params.autocycler_read_type ? params.flye_mode : null) ?: [
        ont_r10: 'nano-hq', ont_r9: 'nano-raw', pacbio_hifi: 'pacbio-hifi', pacbio_clr: 'pacbio-raw'
    ][read_type]
    def flye_mode = mode.startsWith('--') ? mode : "--${mode}"
    def args = params.autocycler_flye_args ?: ''
    def prefix = "${meta.id}.flye.${meta.subsample_id}"
    """
    genome_size=\$(cat ${genome_size_file})
    genome_size_arg=""
    if [ -n "\${genome_size}" ]; then
        genome_size_arg="--genome-size \${genome_size}"
    fi

    flye \\
        ${flye_mode} ${subsample_reads} \\
        \${genome_size_arg} \\
        --threads ${task.cpus} \\
        --out-dir flye_out \\
        ${args}

    cp flye_out/assembly.fasta ${prefix}.assembly.fasta
    if [ -s flye_out/assembly_graph.gfa ]; then
        cp flye_out/assembly_graph.gfa ${prefix}.assembly.gfa
    fi
    cp flye_out/flye.log ${prefix}.flye.log
    cp flye_out/assembly_info.txt ${prefix}.assembly_info.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        flye: \$(flye --version)
    END_VERSIONS
    """
}
