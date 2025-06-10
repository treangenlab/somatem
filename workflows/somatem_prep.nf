#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// -------------------------
// Parameters
// -------------------------
params.input_dir   = 'run1'
params.output_dir  = 'results'
params.threads     = 64
params.maxlength   = 30000
params.minq        = 10
params.minlen      = 250
params.host_index  = 'human-t2t-hla'

// Normalize output_dir - create a derived variable instead of reassigning
def normalized_output_dir = params.output_dir.replaceAll('/+$','')

// -------------------------
// Process: RawNanoPlot1
// -------------------------
process RawNanoPlot1 {
    tag "${sample_id}"
    cpus params.threads
    conda 'bioconda::nanoplot'

    publishDir "${normalized_output_dir}", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    path "${sample_id}_nanoplot1"

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    NanoPlot \\
      -t ${task.cpus} \\
      --fastq "${reads}" \\
      --maxlength ${params.maxlength} \\
      --plots dot \\
      -o "${sample_id}_nanoplot1"
    """
}

// -------------------------
// Process: HostCleanup
// -------------------------
process HostCleanup {
    tag "${sample_id}"
    cpus params.threads
    conda 'bioconda::hostile'
    publishDir "${normalized_output_dir}", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("cleaned_${sample_id}/${sample_id}.clean.fastq.gz"), optional: true

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    mkdir -p cleaned_${sample_id}

    hostile clean \\
      --fastq1 "${reads}" \\
      -t ${task.cpus} \\
      --index ${params.host_index} \\
      --aligner minimap2 \\
      -o cleaned_${sample_id}

    if [[ ! -f cleaned_${sample_id}/${sample_id}.clean.fastq.gz ]]; then
        echo "Warning: cleaned_${sample_id}/${sample_id}.clean.fastq.gz not found → skipping"
        touch cleaned_${sample_id}/${sample_id}.clean.fastq.gz
    fi
    """
}

// -------------------------
// Process: Chopper
// -------------------------
process Chopper {
    tag "${sample_id}"
    cpus params.threads
    conda 'bioconda::chopper'

    publishDir "${normalized_output_dir}", mode: 'copy'

    input:
    tuple val(sample_id), path(cleaned_fastq)

    output:
    tuple val(sample_id), path("chopped_${sample_id}/${sample_id}.chopd.fastq.gz")

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    mkdir -p chopped_${sample_id}

    chopper \\
      -q ${params.minq} \\
      -l ${params.minlen} \\
      -t ${task.cpus} \\
      -i "${cleaned_fastq}" \\
    | gzip > "chopped_${sample_id}/${sample_id}.chopd.fastq.gz"

    if [[ ! -s "chopped_${sample_id}/${sample_id}.chopd.fastq.gz" ]]; then
        echo "Warning: chopped_${sample_id}/${sample_id}.chopd.fastq.gz is empty → creating empty file"
        touch "chopped_${sample_id}/${sample_id}.chopd.fastq.gz"
    fi
    """
}

// -------------------------
// Process: NanoPlot2
// -------------------------
process NanoPlot2 {
    tag "${sample_id}"
    cpus params.threads
    conda 'bioconda::nanoplot'

    publishDir "${normalized_output_dir}", mode: 'copy'

    input:
    tuple val(sample_id), path(chopped_fastq)

    output:
    path "${sample_id}_nanoplot2"

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    NanoPlot \\
      -t ${task.cpus} \\
      --fastq "${chopped_fastq}" \\
      --maxlength ${params.maxlength} \\
      --plots dot \\
      -o "${sample_id}_nanoplot2"
    """
}

// -------------------------
// Process: PrintRunStats
// -------------------------
process PrintRunStats {
    executor 'local'
    
    input:
    val ready

    output:
    stdout

    script:
    """
    #!/usr/bin/env bash
    
    echo "=================================================="
    echo "PIPELINE EXECUTION COMPLETED"
    echo "=================================================="
    echo "Completed at: \$(date +'%d-%b-%Y %H:%M:%S')"
    echo "Pipeline: Long-Read Data Preprocessing"
    echo "Output directory: ${normalized_output_dir}"
    echo "=================================================="
    """
}

// -------------------------
// Workflow Definition
// -------------------------
workflow {
    // Create the base output directory if needed
    new File(normalized_output_dir).mkdirs()

    fastq_ch = Channel
        .fromPath("${params.input_dir}/*.fastq")
        .map { file ->
            def name = file.getName()
            def sample_id = name.replaceFirst(/\.fastq$/,'')
            tuple(sample_id, file)
        }

    // Run initial NanoPlot
    raw_nanoplot_ch = RawNanoPlot1(fastq_ch)

    // Host cleanup - now returns tuple with sample_id
    host_clean_ch = HostCleanup(fastq_ch)

    // Filter out empty results and run Chopper
    chopper_out_ch = Chopper(
        host_clean_ch.filter { sample_id, path -> 
            path.size() > 0 
        }
    )

    // Run final NanoPlot
    final_nanoplot_ch = NanoPlot2(
        chopper_out_ch.filter { sample_id, path -> 
            path.size() > 0 
        }
    )

    // Collect all completion signals and print final stats
    all_complete = raw_nanoplot_ch
        .mix(final_nanoplot_ch)
        .collect()
        .map { "Pipeline completed successfully" }

    PrintRunStats(all_complete) | view
}

// -------------------------
// Workflow Completion Handler
// -------------------------
workflow.onComplete {
    println ""
    println "=================================================="
    println "NEXTFLOW EXECUTION SUMMARY"
    println "=================================================="
    println "Completed at: ${workflow.complete}"
    println "Duration    : ${workflow.duration}"
    println "CPU hours   : ${workflow.stats.computeTimeFmt ?: 'N/A'}"
    println "Succeeded   : ${workflow.stats.succeedCount}"
    if (workflow.errorReport) {
        println "Error report: ${workflow.errorReport}"
    }
    println "=================================================="
    println ""
}