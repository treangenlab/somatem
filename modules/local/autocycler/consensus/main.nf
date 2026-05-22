process AUTOCYCLER_CONSENSUS {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"

    publishDir { "${params.outdir}/isolate_analysis/${meta.id}/autocycler/autocycler_out" },
        mode: params.publish_dir_mode,
        pattern: "autocycler_out/**"
    publishDir { "${params.outdir}/isolate_analysis/${meta.id}/assembly" },
        mode: params.publish_dir_mode,
        pattern: "*.autocycler.{fasta,gfa}"

    input:
    tuple val(meta), path(assemblies)

    output:
    tuple val(meta), path("*.autocycler.fasta"), emit: assembly
    tuple val(meta), path("*.autocycler.gfa")  , emit: graph
    tuple val(meta), path("autocycler_out")    , emit: outdir
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = meta.id
    """
    mkdir -p assemblies
    for assembly in ${assemblies}; do
        cp "\${assembly}" assemblies/
    done

    autocycler compress \\
        -i assemblies \\
        -a autocycler_out

    autocycler cluster \\
        -a autocycler_out

    shopt -s nullglob
    clusters=(autocycler_out/clustering/qc_pass/cluster_*)
    if [ \${#clusters[@]} -eq 0 ]; then
        echo "No Autocycler QC-pass clusters were produced for ${prefix}" >&2
        exit 1
    fi

    for cluster in "\${clusters[@]}"; do
        autocycler trim -c "\${cluster}"
        autocycler resolve -c "\${cluster}"
    done

    autocycler combine \\
        -a autocycler_out \\
        -i autocycler_out/clustering/qc_pass/cluster_*/5_final.gfa

    cp autocycler_out/consensus_assembly.fasta ${prefix}.autocycler.fasta
    cp autocycler_out/consensus_assembly.gfa ${prefix}.autocycler.gfa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        autocycler: \$(autocycler --version 2>&1 | sed 's/^autocycler //')
    END_VERSIONS
    """
}
