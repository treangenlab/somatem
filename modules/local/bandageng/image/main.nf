process BANDAGENG_IMAGE {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3e/3eabbd074e3bc45e2643783450330cae3afc6697fefc635755ab964dc43665a1/data' :
        'community.wave.seqera.io/library/bandage:0.9.0--4f0567049a14ea6d' }"

    input:
    tuple val(meta), path(gfa, arity: '1')
    tuple val(_meta), path(colour_csv, arity: '1') // metadata for node colors

    output:
    tuple val(meta), path('*.png'), emit: png
    tuple val(meta), path('*.svg'), emit: svg
    tuple val("${task.process}"), val('bandageNG'), eval('export QT_QPA_PLATFORM=offscreen; BandageNG --version 2>&1| grep Version | sed "s/^Version: //"'), emit: versions_bandageNG, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def gfa_input = gfa.toString().endsWith('.gz') ? gfa.toString() - ~/\.gz$/ : gfa.toString() // trim the .gz extension
    def decompress = gfa.toString().endsWith('.gz') ? "zcat ${gfa} > ${gfa_input}" : "" // decompress if needed
    def cleanup = gfa.toString().endsWith('.gz') ? "rm ${gfa_input}" : "" // clean up if decompressed

    def colour_by_column = "log_fold_change_0_colour" // choose the column to color by
    def csv_temp="${prefix}_for_bandage_colour.csv"
    """
    # Subset the csv file to 2 columns: col1: node name, col2: ${colour_by_column}
    # If column name "colour" is absent, change "" to "colour"
    csvcut -c 1,"${colour_by_column}" "${colour_csv}" | sed '1s/${colour_by_column}/colour/' > ${csv_temp}

    ${decompress}

    BandageNG image --color ${csv_temp} ${gfa_input} ${prefix}.png $args
    BandageNG image --color ${csv_temp} ${gfa_input} ${prefix}.svg $args

    ${cleanup}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.png
    touch ${prefix}.svg
    """
}
