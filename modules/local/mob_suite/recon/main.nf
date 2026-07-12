process MOB_SUITE_RECON {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"

    publishDir { "${params.outdir}/isolate_analysis/${meta.id}/plasmids" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("${prefix}.contig_report.txt"), emit: contig_report
    tuple val(meta), path("${prefix}.mobtyper_results.txt"), optional: true, emit: mobtyper
    tuple val(meta), path("${prefix}.mge_report.txt"), optional: true, emit: mge_report
    tuple val(meta), path("${prefix}.plasmids/*.fasta"), optional: true, emit: plasmids
    path "versions.yml", emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def database_dir = params.mob_suite_db ?: "${params.db_base_dir}/mob_suite"
    """
    mkdir -p '${database_dir}' ${prefix}.plasmids

    if [[ ! -s '${database_dir}/clusters.txt' ]]; then
        echo 'Initializing the MOB-suite plasmid databases. This is required only once.' >&2
        mob_init --database_directory '${database_dir}'
    fi

    mob_recon \\
        --infile ${assembly} \\
        --outdir mob_recon_out \\
        --database_directory '${database_dir}' \\
        --num_threads ${task.cpus}

    cp mob_recon_out/contig_report.txt ${prefix}.contig_report.txt
    if [[ -f mob_recon_out/mobtyper_results.txt ]]; then
        cp mob_recon_out/mobtyper_results.txt ${prefix}.mobtyper_results.txt
    elif [[ -f mob_recon_out/mobtyper_results ]]; then
        cp mob_recon_out/mobtyper_results ${prefix}.mobtyper_results.txt
    fi
    if [[ -f mob_recon_out/mge.report.txt ]]; then
        cp mob_recon_out/mge.report.txt ${prefix}.mge_report.txt
    fi
    shopt -s nullglob
    for plasmid in mob_recon_out/plasmid_*.fasta; do
        cp "\${plasmid}" ${prefix}.plasmids/
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mob_suite: \$(mob_recon --version 2>&1 | head -1)
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}.plasmids
    printf 'sample_id\tcontig_id\tmolecule_type\tsize\tcircularity_status\tpredicted_mobility\n${meta.id}\tcontig_1\tChromosome\t1\tnot tested\t-\n' > ${prefix}.contig_report.txt
    touch ${prefix}.mobtyper_results.txt ${prefix}.mge_report.txt
    echo '>${meta.id}_plasmid' > ${prefix}.plasmids/plasmid_1.fasta
    echo 'N' >> ${prefix}.plasmids/plasmid_1.fasta
    printf '"${task.process}":\n    mob_suite: "stub"\n' > versions.yml
    """
}
