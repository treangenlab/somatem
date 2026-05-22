process SOMATEM_SUMMARY_REPORT {
    tag "$workflow_slug"
    label 'process_single'

    conda "conda-forge::python=3.12"

    publishDir { "${params.outdir}/summary_reports" }, mode: params.publish_dir_mode, pattern: "*.html"

    input:
    val workflow_slug
    val workflow_label
    path samplesheet
    path report_files, stageAs: "report_inputs/file??/*"

    output:
    path "${workflow_slug}_summary_report.html", emit: html
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def pipeline_version = workflow.manifest.version ?: "unknown"
    def run_name = workflow.runName ?: "unknown"
    def workflow_slug_q = "'" + (workflow_slug == null ? "" : workflow_slug.toString()).replace("'", "'\"'\"'") + "'"
    def workflow_label_q = "'" + (workflow_label == null ? "" : workflow_label.toString()).replace("'", "'\"'\"'") + "'"
    def analysis_type_q = "'" + (params.analysis_type == null ? "" : params.analysis_type.toString()).replace("'", "'\"'\"'") + "'"
    def data_type_q = "'" + (params.data_type == null ? "" : params.data_type.toString()).replace("'", "'\"'\"'") + "'"
    def sample_environment_q = "'" + (params.sample_environment == null ? "" : params.sample_environment.toString()).replace("'", "'\"'\"'") + "'"
    def sequencing_technology_q = "'" + (params.sequencing_technology == null ? "" : params.sequencing_technology.toString()).replace("'", "'\"'\"'") + "'"
    def pipeline_version_q = "'" + (pipeline_version == null ? "" : pipeline_version.toString()).replace("'", "'\"'\"'") + "'"
    def run_name_q = "'" + (run_name == null ? "" : run_name.toString()).replace("'", "'\"'\"'") + "'"
    """
    python3 ${projectDir}/bin/somatem_report.py \\
        --workflow ${workflow_slug_q} \\
        --workflow-label ${workflow_label_q} \\
        --samplesheet ${samplesheet} \\
        --input-dir report_inputs \\
        --output ${workflow_slug}_summary_report.html \\
        --analysis-type ${analysis_type_q} \\
        --data-type ${data_type_q} \\
        --sample-environment ${sample_environment_q} \\
        --sequencing-technology ${sequencing_technology_q} \\
        --pipeline-version ${pipeline_version_q} \\
        --run-name ${run_name_q}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        somatem_report: "1.0.0"
        python: \$(python3 --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    cat <<-END_REPORT > ${workflow_slug}_summary_report.html
    <!doctype html><html><body><h1>${workflow_label}</h1><p>Stub somatem summary report.</p></body></html>
    END_REPORT

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        somatem_report: "1.0.0"
        python: "stub"
    END_VERSIONS
    """
}
