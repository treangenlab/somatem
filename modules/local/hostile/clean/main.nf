process HOSTILE_CLEAN {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/hostile:2.0.2--a7f5e5d341b6b94b"

    input:
    tuple val(meta), path(reads, stageAs: 'reads.fastq.gz', arity: '1')
    path index, stageAs: 'host.mmi'

    output:
    tuple val(meta), path('*.hostile.clean.fastq.gz'), emit: fastq
    tuple val(meta), path('*.hostile.summary.json'), emit: json
    path 'versions.yml', emit: versions

    script:
    def sample = meta.id.replace("'", "'\\''")
    """
    export HOSTILE_CACHE_DIR=./hostile-cache
    hostile clean --aligner minimap2 --index host.mmi --fastq1 reads.fastq.gz \\
        --threads ${task.cpus} --airplane --output . > '${sample}.hostile.summary.json'
    mv reads.clean.fastq.gz '${sample}.hostile.clean.fastq.gz'

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hostile: \$(hostile --version)
        minimap2: \$(minimap2 --version)
    END_VERSIONS
    """

    stub:
    def sample = meta.id.replace("'", "'\\''")
    """
    gzip -c </dev/null > '${sample}.hostile.clean.fastq.gz'
    echo '[]' > '${sample}.hostile.summary.json'
    echo 'hostile: stub' > versions.yml
    """
}
