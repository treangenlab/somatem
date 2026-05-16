process PIGEON_ITERATIVE_BINNING {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(gfa), path(assembly), path(bam)

    output:
    tuple val(meta), path("iterative_binning/final_bins/*"), emit: final_bins
    tuple val(meta), path("iterative_binning/candidate_manifest.tsv"), emit: candidate_manifest
    tuple val(meta), path("iterative_binning/selection/selected_iterations.tsv"), emit: selected_tsv
    tuple val(meta), path("iterative_binning/selection/selection_trajectory.tsv"), emit: trajectory_tsv
    tuple val(meta), path("iterative_binning/selection/selection_summary.json"), emit: selection_summary
    tuple val(meta), path("iterative_binning/iterative_summary.json"), emit: iterative_summary
    tuple val(meta), path("iterative_binning/command_log.tsv"), emit: command_log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    pigeon_iterate.py \\
        --gfa ${gfa} \\
        --assembly ${assembly} \\
        --bam ${bam} \\
        --outdir iterative_binning \\
        --sample-id "${meta.id}" \\
        --threads $task.cpus \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pigeon_iterate: "1.0.0"
        python: \$(python --version 2>&1 | sed 's/Python //')
        semibin: \$(SemiBin2 --version 2>&1 || true)
        metabat2: \$(metabat2 --version 2>&1 | head -n 1 || true)
        vamb: \$(vamb --version 2>&1 | head -n 1 || true)
        das_tool: \$(DAS_Tool --version 2>&1 | head -n 1 || true)
    END_VERSIONS
    """

    stub:
    """
    mkdir -p iterative_binning/final_bins
    mkdir -p iterative_binning/selection
    touch iterative_binning/final_bins/stub.fa
    cat <<-END_MANIFEST > iterative_binning/candidate_manifest.tsv
    sample	binner	iteration	seed	target_fasta	residual_fasta	metrics	bins
    ${meta.id}	semibin2	1	42	assembly.fa	residual.fa	metrics.json	bins
    END_MANIFEST
    cat <<-END_SELECTED > iterative_binning/selection/selected_iterations.tsv
    sample	binner	iteration	seed	selection_score	decision	pam	bin_explained_fraction	both_fraction	unexplained_fraction	residual_contig_fraction	residual_kmer_fraction	trajectory_stop_iteration	trajectory_stop_reason	target_fasta	residual_fasta	metrics	bins	selected_bins_dir
    ${meta.id}	semibin2	1	42	1.000000	stub	1.0	1.0	1.0	0.0	0.0	0.0	1	stub	assembly.fa	residual.fa	metrics.json	bins	selected_bins
    END_SELECTED
    cat <<-END_TRAJECTORY > iterative_binning/selection/selection_trajectory.tsv
    sample	binner	iteration	seed	score	delta_score	best_score_so_far	no_gain_count	gain_state	decision	pam	bin_explained_fraction	both_fraction	unexplained_fraction	residual_contig_fraction	residual_kmer_fraction	target_fasta	residual_fasta
    ${meta.id}	semibin2	1	42	1.0		1.0	0	improved	stub	1.0	1.0	1.0	0.0	0.0	0.0	assembly.fa	residual.fa
    END_TRAJECTORY
    cat <<-END_SELECTION_SUMMARY > iterative_binning/selection/selection_summary.json
    {"schema_version":"pigeon.selection.v1","score_mode":"stub","n_candidates":1,"n_selected":1,"selected":[]}
    END_SELECTION_SUMMARY
    cat <<-END_ITERATIVE_SUMMARY > iterative_binning/iterative_summary.json
    {"schema_version":"pigeon.iterative_binning.v1","sample_id":"${meta.id}","n_candidates":1}
    END_ITERATIVE_SUMMARY
    cat <<-END_COMMAND_LOG > iterative_binning/command_log.tsv
    cwd	command
    stub	stub
    END_COMMAND_LOG

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pigeon_iterate: "1.0.0"
        python: "stub"
    END_VERSIONS
    """
}
