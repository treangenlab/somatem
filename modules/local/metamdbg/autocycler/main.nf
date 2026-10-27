process METAMDBG_AUTOCYCLER {
    tag "${meta.id}:${meta.subsample_id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"

    publishDir { "${params.outdir}/isolate_analysis/${meta.id}/autocycler/assemblies" },
        mode: params.publish_dir_mode,
        pattern: "*.assembly.{fasta,gfa,log}"

    input:
    tuple val(meta), path(subsample_reads), path(genome_size_file), val(species), val(read_type)

    output:
    tuple val(meta), path("*.assembly.fasta"), emit: assembly
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = "${meta.id}.metamdbg.${meta.subsample_id}.assembly"
    def args = params.autocycler_metamdbg_args ?: ''
    def extra_args = args ? "--args ${args}" : ''
    """
    autocycler helper metamdbg \\
        --reads "${subsample_reads}" \\
        --out_prefix "${prefix}" \\
        --threads ${task.cpus} \\
        --read_type ${read_type} \\
        --genome_size "\$(cat '${genome_size_file}')" \\
        --dir metamdbg_out \\
        ${extra_args}

    test -s "${prefix}.fasta"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        autocycler: \$(autocycler --version 2>&1 | sed 's/^autocycler //')
        metamdbg: \$(metaMDBG 2>&1 | sed -n 's/.*Version: //p')
    END_VERSIONS
    """
}
