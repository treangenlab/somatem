process SEQSCREEN_REFERENCE_INFERENCE {
    tag "$meta.id"
    label 'process_medium'

    conda "${projectDir}/modules/local/seqscreen/reference_inference/environment.yml"

    input:
    tuple val(meta), path(fasta), path(report), path(taxonomy)
    path db
    path assets

    output:
    tuple val(meta), path("reference_inference"), emit: outdir
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def online = params.seqscreen_online ? '--online' : ''
    """
    mkdir -p reference_inference/taxonomic_identification/taxonomic_assignment inference_working
    cp ${taxonomy} reference_inference/taxonomic_identification/taxonomic_assignment/taxonomic_results.txt

    python3 ${assets}/scripts/reference_inference_short_read.py \\
        --fasta1=${fasta} \\
        --output=reference_inference \\
        --working=inference_working \\
        --databases=${db} \\
        --threads=${task.cpus} \\
        ${online}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqscreen: "4.5"
    END_VERSIONS
    """

    stub:
    """
    mkdir -p reference_inference
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqscreen: "stub"
    END_VERSIONS
    """
}
