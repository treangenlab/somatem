process AGB {
    tag "${meta.id}"
    label 'process_low'

    publishDir { "${params.outdir}/assembly_mags/${meta.id}/assembly_graph" },
        mode: params.publish_dir_mode,
        pattern: 'agb_output/**'

    conda 'almiheenko::agb'

    input:
    tuple val(meta), path(gfa), path(graph), path(info)

    output:
    tuple val(meta), path('agb_output/viewer.html'), emit: assembly_graph
    path 'versions.yml', emit: versions

    script:
    """
    python ${projectDir}/bin/run_agb.py \
        --gfa ${gfa} \
        --graph ${graph} \
        --info ${info} \
        --assembler flye \
        --threads ${task.cpus}
    """

    stub:
    """
    mkdir -p agb_output
    touch agb_output/viewer.html
    echo 'AGB:' > versions.yml
    echo '    agb: "stub"' >> versions.yml
    """
}
