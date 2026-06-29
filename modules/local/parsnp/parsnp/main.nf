process PARSNP {
    tag "parsnp"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/parsnp:2.1.4--pyhdfd78af_0':
        'quay.io/biocontainers/parsnp:2.1.4--pyhdfd78af_0' }"

    input:
    path reference
    path genomes

    output:
    path "parsnp_out"                  , emit: outdir
    path "parsnp_out/parsnp.ggr"       , optional: true, emit: gingr_archive
    path "parsnp_out/parsnp.xmfa"      , optional: true, emit: xmfa
    path "parsnp_out/parsnp.tree"      , optional: true, emit: tree
    path "parsnp_out/parsnp.snps.mblocks", emit: core_snps_fasta
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mkdir -p genomes parsnp_out

    for genome in ${genomes.join(' ')}; do
        ln -sf ../\${genome} genomes/\$(basename \${genome})
    done

    parsnp \\
        -r ${reference} \\
        -d genomes \\
        -o parsnp_out \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        parsnp: \$(parsnp -V 2>&1 | head -n 1 || true)
    END_VERSIONS
    """

    stub:
    """
    mkdir -p parsnp_out
    cat <<-END_FASTA > parsnp_out/parsnp.snps.mblocks
    >reference
    A
    >sample
    G
    END_FASTA
    echo "(reference,sample);" > parsnp_out/parsnp.tree
    touch parsnp_out/parsnp.ggr
    touch parsnp_out/parsnp.xmfa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        parsnp: "stub"
    END_VERSIONS
    """
}
