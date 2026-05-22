process AUTOCYCLER_SUBSAMPLE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"

    publishDir { "${params.outdir}/isolate_analysis/${meta.id}/autocycler/subsamples" },
        mode: params.publish_dir_mode,
        pattern: "subsampled_reads/*.fastq"

    input:
    tuple val(meta), path(long_reads), val(genome_size), val(species)

    output:
    tuple val(meta), path("subsampled_reads/sample_*.fastq"), path("genome_size.txt"), val(species), emit: reads
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def count = params.autocycler_subsamples ?: 4
    def min_read_depth = params.autocycler_min_read_depth ?: params.autocycler_subsample_depth ?: 25
    def seed = params.autocycler_seed ?: 0
    def genome_size_string = genome_size ?: ''
    """
    genome_size="${genome_size_string}"
    if [ -z "\${genome_size}" ] || [ "\${genome_size}" = "null" ]; then
        genome_size=\$(autocycler helper genome_size --reads ${long_reads} --threads ${task.cpus})
    fi
    printf "%s\\n" "\${genome_size}" > genome_size.txt

    autocycler subsample \\
        --reads ${long_reads} \\
        --out_dir subsampled_reads \\
        --genome_size "\${genome_size}" \\
        --count ${count} \\
        --min_read_depth ${min_read_depth} \\
        --seed ${seed}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        autocycler: \$(autocycler --version 2>&1 | sed 's/^autocycler //')
    END_VERSIONS
    """
}
