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

// Normalize output_dir exactly once
params.output_dir = params.output_dir.replaceAll('/+$','')

// -------------------------
// Process: RawNanoPlot1 (fixed)
// -------------------------
process RawNanoPlot1 {
    tag "${sample_id}"
    cpus params.threads
    conda 'bioconda::nanoplot'

    // publishDir must use named arguments
    publishDir path: "${params.output_dir}", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    // Use ** instead of **/* to capture files directly inside the folder
    path "${sample_id}_nanoplot1/**"

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
// Process: HostCleanup (unchanged)
// -------------------------
process HostCleanup {
    tag "${sample_id}"
    cpus params.threads

    conda """
      channels:
        - bioconda
        - conda-forge
        - defaults
      dependencies:
        - hostile
        - minimap2
    """

    publishDir path: "${params.output_dir}", mode: 'copy', pattern: "cleaned_${sample_id}/*"

    input:
    tuple val(sample_id), path(reads)

    output:
    path "cleaned_${sample_id}/${sample_id}.clean.fastq.gz", optional: true

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
        exit 0
    fi
    """
}

// -------------------------
// Process: Chopper (unchanged)
// -------------------------
process Chopper {
    tag "${sample_id}"
    cpus params.threads
    conda 'bioconda::chopper'

    publishDir path: "${params.output_dir}", mode: 'copy', pattern: "chopped_${sample_id}/*"

    input:
    tuple val(sample_id), path(cleaned_fastq)

    output:
    path "chopped_${sample_id}/${sample_id}.chopd.fastq.gz"

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
        echo "Warning: chopped_${sample_id}/${sample_id}.chopd.fastq.gz is empty → skipping"
        exit 0
    fi
    """
}

// -------------------------
// Process: NanoPlot2 (example, also fixed if needed)
// -------------------------
process NanoPlot2 {
    tag "${sample_id}"
    cpus params.threads
    conda 'bioconda::nanoplot'

    publishDir path: "${params.output_dir}", mode: 'copy'

    input:
    tuple val(sample_id), path(chopped_fastq)

    output:
    path "${sample_id}_nanoplot2/**"

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
// Workflow Definition
// -------------------------
workflow {
    // Create the base output directory if needed
    new File(params.output_dir).mkdirs()

    fastq_ch = Channel
        .fromPath("${params.input_dir}/*.fastq")
        .map { file ->
            def name = file.getName()                 // e.g. "45_2.fastq"
            def sample_id = name.replaceFirst(/\.fastq$/,'')
            tuple(sample_id, file)
        }

    RawNanoPlot1(fastq_ch)

    host_clean_ch = HostCleanup(fastq_ch)

    chopper_input_ch = host_clean_ch
        .map { path ->
            def fname = path.getName()                // e.g. "45_2.clean.fastq.gz"
            def sample_id = fname.replaceFirst(/\.clean\.fastq\.gz$/,'')
            tuple(sample_id, path)
        }

    chopper_out_ch = Chopper(chopper_input_ch)

    nanoplot2_input_ch = chopper_out_ch
        .map { path ->
            def fname = path.getName()                // e.g. "45_2.chopd.fastq.gz"
            def sample_id = fname.replaceFirst(/\.chopd\.fastq\.gz$/,'')
            tuple(sample_id, path)
        }

    NanoPlot2(nanoplot2_input_ch)
}

