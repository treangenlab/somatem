#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/* SeqScreen that composes the repo's .nf workflows:

- Initialization
- SeqMapper
- Taxonomic Identification (Fast or Sensitive mode based on cfg.mode parameter)
- Functional Annotation
- Report Generation

Usage:
nextflow run seqscreen.nf \
    --input_dir /path/to/input/directory \
    --databases /path/to/seqscreen/databases \
    --bin_dir /path/to/seqscreen/bin \
    --workflows_dir /path/to/seqscreen/workflows \
    --module_dir /path/to/seqscreen/modules
*/

// Help message
def helpMessage() {
    log.info"""
Usage:
    nextflow run seqscreen.nf [options]

Required Arguments:
    --input_dir     Path to input directory containing FASTA/FASTQ files
    --databases     Path to SeqScreen database directory
    --bin_dir       Path to SeqScreen bin directory
    --workflows_dir Path to SeqScreen workflows directory  
    --module_dir    Path to SeqScreen modules directory

Optional Arguments:
    --working       Working directory (default: work)
    --mode          Analysis mode: fast or sensitive (default: fast)
    --threads       Number of threads (default: 1)
    --evalue        E-value threshold (default: 10)
    --slurm         Enable SLURM execution (default: false)
    --version       SeqScreen version (default: 1.0)
    --prefix        Output file prefix (default: "")
    --max_rapsearch_threads  Maximum threads for rapsearch2 (default: 16)
    --min_file_size_mb       Minimum file size in MB for full resources (default: 1)

Supported File Formats:
    FASTA: .fa, .fasta, .fas, .fna (and their .gz versions)
    FASTQ: .fq, .fastq (and their .gz versions)

Profiles:
    -profile standard    Run locally (default)
    -profile slurm       Run on SLURM cluster

Example:
    nextflow run seqscreen.nf \\
        --input_dir /path/to/sequences/ \\
        --databases /path/to/db \\
        --workflows_dir /path/to/workflows \\
        --bin_dir /path/to/bin \\
        --module_dir /path/to/modules \\
        --threads 8
""".stripIndent()
}

// Show help message if requested
if (params.help) {
    helpMessage()
    exit 0
}

//
// Process to convert FASTQ files to FASTA format using seqkit fq2fa
//
process SEQKIT_FQ2FA {
    tag "${meta.id}"
    label 'process_low'
    conda 'bioconda::seqkit'

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("*.fasta"), emit: fasta

    script:
    def output_name = "${meta.id}.fasta"
    def is_gzipped = fastq.name.endsWith('.gz')
    
    """
    if [ "${is_gzipped}" = "true" ]; then
        # Handle gzipped FASTQ files
        seqkit fq2fa ${fastq} -o ${output_name}
    else
        # Handle uncompressed FASTQ files
        seqkit fq2fa ${fastq} -o ${output_name}
    fi
    """
}

//
// Wrapper processes that set parameters and call the existing workflows
//

process RUN_INITIALIZE {
    tag "${meta.id}"
    label 'process_low'
    conda 'bioconda::seqscreen'

    input:
    tuple val(meta), path(fasta)
    path db_dir
    val cfg

    output:
    tuple val(meta), path(fasta), emit: initialized

    script:
    def working_dir = cfg?.working ?: "work"
    def log_file = cfg?.log ?: "/dev/null"
    def slurm_flag = cfg?.slurm ? "--slurm" : ""
    def sensitive_flag = cfg?.mode == "sensitive" ? "--sensitive" : ""
    def workflows_dir = cfg?.workflows_dir ?: "."
    def bin_dir = cfg?.bin_dir ?: "."
    def module_dir = cfg?.module_dir ?: "."

    """
    # Set up parameters and run initialize workflow
    export NXF_WORK=\${PWD}/work_init

    # Set environment variables that the sub-subworkflow expects
    export BIN_DIR=${bin_dir}
    export MODULE_DIR=${module_dir}

    nextflow run ${workflows_dir}/initialize.nf \\
        --fasta ${fasta} \\
        --working ${working_dir} \\
        --bin_dir ${bin_dir} \\
        --module_dir ${module_dir} \\
        ${slurm_flag} \\
        ${sensitive_flag} \\
        -work-dir \${NXF_WORK}
    """
}

process RUN_SEQMAPPER {
    tag "${meta.id}"
    label 'process_medium'
    conda 'bioconda::seqscreen'

    input:
    tuple val(meta), path(fasta)
    path db_dir
    val cfg

    output:
    tuple val(meta), path(fasta), emit: mapped

    script:
    def working_dir = cfg?.working ?: "work"
    def log_file = cfg?.log ?: "/dev/null"
    def threads = cfg?.threads ?: 1
    def slurm_flag = cfg?.slurm ? "--slurm" : ""
    def hmmscan_flag = cfg?.hmmscan ? "--hmmscan" : ""
    def workflows_dir = cfg?.workflows_dir ?: "."
    def bin_dir = cfg?.bin_dir ?: "."
    def module_dir = cfg?.module_dir ?: "."
    def max_rapsearch_threads = cfg?.max_rapsearch_threads ?: 16
    def min_file_size_mb = cfg?.min_file_size_mb ?: 1

    """
    # Set up parameters and run seqmapper workflow
    export NXF_WORK=\${PWD}/work_seqmapper

    # Set environment variables that the sub-subworkflow expects
    export BIN_DIR=${bin_dir}
    export MODULE_DIR=${module_dir}

    # Log file size for debugging
    FILE_SIZE=\$(stat -c%s "${fasta}")
    FILE_SIZE_MB=\$((\$FILE_SIZE / 1024 / 1024))
    echo "SeqMapper input file size: \$FILE_SIZE bytes (\$FILE_SIZE_MB MB)" | tee -a ${log_file}

    nextflow run ${workflows_dir}/seqmapper.nf \\
        --fasta ${fasta} \\
        --working ${working_dir} \\
        --databases ${db_dir} \\
        --bin_dir ${bin_dir} \\
        --module_dir ${module_dir} \\
        --threads ${threads} \\
        --max_rapsearch_threads ${max_rapsearch_threads} \\
        --min_file_size_mb ${min_file_size_mb} \\
        ${slurm_flag} \\
        ${hmmscan_flag} \\
        -work-dir \${NXF_WORK}
    """
}

process RUN_TAXONOMIC_IDENTIFICATION_FAST {
    tag "${meta.id}"
    label 'process_high'
    conda 'bioconda::seqscreen'

    input:
    tuple val(meta), path(fasta)
    path db_dir
    val cfg

    output:
    tuple val(meta), path(fasta), emit: taxonomic_results

    script:
    def working_dir = cfg?.working ?: "work"
    def log_file = cfg?.log ?: "/dev/null"
    def threads = cfg?.threads ?: 1
    def evalue = cfg?.evalue ?: 10
    def taxlimit = cfg?.taxlimit ?: 25
    def splitby = cfg?.splitby ?: 100000
    def slurm_flag = cfg?.slurm ? "--slurm" : ""
    def workflows_dir = cfg?.workflows_dir ?: "."
    def bin_dir = cfg?.bin_dir ?: "."
    def module_dir = cfg?.module_dir ?: "."

    """
    # Set up parameters and run fast taxonomic identification workflow
    export NXF_WORK=\${PWD}/work_tax_fast

    # Set environment variables that the sub-subworkflow expects
    export BIN_DIR=${bin_dir}
    export MODULE_DIR=${module_dir}

    nextflow run ${workflows_dir}/taxonomic_identification_fast.nf \\
        --fasta ${fasta} \\
        --working ${working_dir} \\
        --databases ${db_dir} \\
        --bin_dir ${bin_dir} \\
        --module_dir ${module_dir} \\
        --threads ${threads} \\
        --evalue ${evalue} \\
        --taxlimit ${taxlimit} \\
        --splitby ${splitby} \\
        ${slurm_flag} \\
        -work-dir \${NXF_WORK}
    """
}

process RUN_TAXONOMIC_IDENTIFICATION_SENSITIVE {
    tag "${meta.id}"
    label 'process_high'
    conda 'bioconda::seqscreen'

    input:
    tuple val(meta), path(fasta)
    path db_dir
    val cfg

    output:
    tuple val(meta), path(fasta), emit: taxonomic_results

    script:
    def working_dir = cfg?.working ?: "work"
    def log_file = cfg?.log ?: "/dev/null"
    def threads = cfg?.threads ?: 1
    def evalue = cfg?.evalue ?: 10
    def slurm_flag = cfg?.slurm ? "--slurm" : ""
    def blastn_flag = cfg?.blastn ? "--blastn" : ""
    def workflows_dir = cfg?.workflows_dir ?: "."
    def bin_dir = cfg?.bin_dir ?: "."
    def module_dir = cfg?.module_dir ?: "."

    """
    # Set up parameters and run sensitive taxonomic identification workflow
    export NXF_WORK=\${PWD}/work_tax_sensitive

    # Set environment variables that the sub-subworkflow expects
    export BIN_DIR=${bin_dir}
    export MODULE_DIR=${module_dir}

    nextflow run ${workflows_dir}/taxonomic_identification_sensitive.nf \\
        --fasta ${fasta} \\
        --working ${working_dir} \\
        --databases ${db_dir} \\
        --bin_dir ${bin_dir} \\
        --module_dir ${module_dir} \\
        --threads ${threads} \\
        --evalue ${evalue} \\
        ${slurm_flag} \\
        ${blastn_flag} \\
        -work-dir \${NXF_WORK}
    """
}

process RUN_FUNCTIONAL_ANNOTATION {
    tag "${meta.id}"
    label 'process_high'
    conda 'bioconda::seqscreen'

    input:
    tuple val(meta), path(fasta)  // This comes from taxonomic identification now
    path db_dir
    val cfg

    output:
    tuple val(meta), path(fasta), emit: functional_results

    script:
    def working_dir = cfg?.working ?: "work"
    def log_file = cfg?.log ?: "/dev/null"
    def threads = cfg?.threads ?: 1
    def evalue = cfg?.evalue ?: (cfg?.mode == "sensitive" ? 10 : 30)
    def slurm_flag = cfg?.slurm ? "--slurm" : ""
    def ancestral_flag = cfg?.ancestral ? "--ancestral" : ""
    def sensitive_flag = cfg?.mode == "sensitive" ? "--sensitive" : ""
    def includecc_flag = cfg?.includecc ? "--includecc" : ""
    def bitscore_cutoff = cfg?.bitscore_cutoff ?: 5
    def workflows_dir = cfg?.workflows_dir ?: "."
    def bin_dir = cfg?.bin_dir ?: "."
    def module_dir = cfg?.module_dir ?: "."

    """
    # Set up parameters and run functional annotation workflow
    export NXF_WORK=\${PWD}/work_functional

    # Set environment variables that the sub-subworkflow expects
    export BIN_DIR=${bin_dir}
    export MODULE_DIR=${module_dir}

    nextflow run ${workflows_dir}/functional_annotation.nf \\
        --fasta ${fasta} \\
        --working ${working_dir} \\
        --databases ${db_dir} \\
        --bin_dir ${bin_dir} \\
        --module_dir ${module_dir} \\
        --threads ${threads} \\
        --evalue ${evalue} \\
        --bitscore_cutoff ${bitscore_cutoff} \\
        ${slurm_flag} \\
        ${ancestral_flag} \\
        ${sensitive_flag} \\
        ${includecc_flag} \\
        -work-dir \${NXF_WORK}
    """
}

process RUN_REPORT_GENERATION {
    tag "${meta.id}"
    label 'process_low'
    publishDir "${params.working}/results", mode: 'copy'
    conda 'bioconda::seqscreen'

    input:
    tuple val(meta), path(fasta)  // From taxonomic identification
    val functional_complete       // Just a signal that functional annotation is done
    path db_dir
    val cfg

    output:
    tuple val(meta), path("report_output/*"), emit: html
    tuple val(meta), path("report_output/*.tsv"), emit: tables
    tuple val(meta), path("report_output/*.json"), emit: json

    script:
    def working_dir = cfg?.working ?: "work"
    def log_file = cfg?.log ?: "/dev/null"
    def version = cfg?.version ?: "1.0"
    def prefix = cfg?.prefix ? "${cfg.prefix}_" : ""
    def slurm_flag = cfg?.slurm ? "--slurm" : ""
    def sensitive_flag = cfg?.mode == "sensitive" ? "--sensitive" : ""
    def ancestral_flag = cfg?.ancestral ? "--ancestral" : ""
    def hmmscan_flag = cfg?.hmmscan ? "--hmmscan" : ""
    def blastn_flag = cfg?.blastn ? "--blastn" : ""
    def includecc_flag = cfg?.includecc ? "--includecc" : ""
    def evalue_flag = cfg?.evalue ? "--evalue ${cfg.evalue}" : ""
    def bitscore_flag = cfg?.bitscore_cutoff ? "--bitscore ${cfg.bitscore_cutoff}" : ""
    def taxlimit_flag = cfg?.taxlimit ? "--taxlimit ${cfg.taxlimit}" : ""
    def splitby_flag = cfg?.splitby ? "--splitby ${cfg.splitby}" : ""
    def workflows_dir = cfg?.workflows_dir ?: "."
    def bin_dir = cfg?.bin_dir ?: "."
    def module_dir = cfg?.module_dir ?: "."

    """
    # Set up parameters and run report generation workflow
    export NXF_WORK=\${PWD}/work_report

    # Set environment variables that the sub-subworkflow expects
    export BIN_DIR=${bin_dir}
    export MODULE_DIR=${module_dir}

    nextflow run ${workflows_dir}/report_generation.nf \\
        --fasta ${fasta} \\
        --working ${working_dir} \\
        --databases ${db_dir} \\
        --bin_dir ${bin_dir} \\
        --module_dir ${module_dir} \\
        --version ${version} \\
        --prefix ${prefix} \\
        ${slurm_flag} \\
        ${sensitive_flag} \\
        ${ancestral_flag} \\
        ${hmmscan_flag} \\
        ${blastn_flag} \\
        ${includecc_flag} \\
        ${evalue_flag} \\
        ${bitscore_flag} \\
        ${taxlimit_flag} \\
        ${splitby_flag} \\
        -work-dir \${NXF_WORK}

    # Create output directory and copy results
    mkdir -p report_output

    # Copy HTML reports if they exist
    if [ -d "${working_dir}/report_generation/${prefix}seqscreen_html_report" ]; then
        cp -r ${working_dir}/report_generation/${prefix}seqscreen_html_report/* report_output/
    fi

    # Copy TSV reports if they exist
    if [ -f "${working_dir}/report_generation/${prefix}seqscreen_report.tsv" ]; then
        cp ${working_dir}/report_generation/${prefix}seqscreen_report.tsv report_output/
    fi

    # Create a dummy JSON file if none exists (since we emit json but the original only creates TSV)
    if [ ! -f "report_output/*.json" ]; then
        echo '{"status": "completed", "sample": "'${meta.id}'"}' > report_output/summary.json
    fi
    """
}

workflow SEQSCREEN {

    take:
    seqs_ch
    db_dir
    cfg

    main:
    //
    // Adapter layer
    // Make tiny, explicit tuples so the ports line up regardless of your upstream.
    //
    ch_input = seqs_ch
    ch_db    = db_dir
    ch_cfg   = cfg

    //
    // 1) Initialization
    // Initialize the SeqScreen workflow and validate input
    //
    RUN_INITIALIZE(ch_input, ch_db, ch_cfg)
    init_seqs = RUN_INITIALIZE.out.initialized

    //
    // 2) Mapping
    // Run sequence mapping against databases
    //
    RUN_SEQMAPPER(init_seqs, ch_db, ch_cfg)
    mapped = RUN_SEQMAPPER.out.mapped

    //
    // 3) Taxonomic identification
    // Call different sub-subworkflows based on user parameters (mode: fast/sensitive)
    // IMPORTANT: This must complete before functional annotation starts
    //

    // Determine which taxonomic identification workflow to run based on mode
    if (ch_cfg.mode == "sensitive") {
        // Run sensitive taxonomic identification
        RUN_TAXONOMIC_IDENTIFICATION_SENSITIVE(
            mapped,
            ch_db,
            ch_cfg
        )
        tax_calls = RUN_TAXONOMIC_IDENTIFICATION_SENSITIVE.out.taxonomic_results
    } else {
        // Run fast taxonomic identification (default)
        RUN_TAXONOMIC_IDENTIFICATION_FAST(
            mapped,
            ch_db,
            ch_cfg
        )
        tax_calls = RUN_TAXONOMIC_IDENTIFICATION_FAST.out.taxonomic_results
    }

    //
    // 4) Functional annotation
    // FIXED: Now takes input from taxonomic identification (tax_calls) instead of mapped
    // This ensures taxonomic identification completes before functional annotation starts
    //
    RUN_FUNCTIONAL_ANNOTATION(tax_calls, ch_db, ch_cfg)
    func_calls = RUN_FUNCTIONAL_ANNOTATION.out.functional_results

    //
    // 5) Report generation
    // FIXED: Now takes only tax_calls (with FASTA) and a completion signal from functional annotation
    // This avoids the file name collision while ensuring both processes complete first
    //
    RUN_REPORT_GENERATION(
        tax_calls,                    // Contains meta and fasta
        func_calls.map { it -> true }, // Just a completion signal, no files
        ch_db, 
        ch_cfg
    )

    //
    // 6) Emit canonical ports you can hook into MultiQC/Dash/etc.
    //
    emit:
    report_html_ch = RUN_REPORT_GENERATION.out.html
    tables_ch      = RUN_REPORT_GENERATION.out.tables
    json_ch        = RUN_REPORT_GENERATION.out.json
}

//
// Helper function to extract proper file ID from filename
//
def getFileId(file) {
    def name = file.name
    // Remove common extensions in order: .gz, then .fasta/.fa/.fastq/.fq/.fas/.fna
    name = name.replaceAll(/\.gz$/, '')
    name = name.replaceAll(/\.(fasta|fa|fastq|fq|fas|fna)$/, '')
    return name
}

//
// Main workflow - Entry point for Nextflow execution
//
workflow {

    // Validate required parameters
    if (!params.input_dir) {
        error "Please provide an input directory with --input_dir"
    }
    if (!params.databases) {
        error "Please provide a database directory with --databases"
    }
    if (!params.bin_dir) {
        error "Please provide the SeqScreen bin directory with --bin_dir (contains helper scripts like validate_fasta.pl)"
    }
    if (!params.workflows_dir) {
        error "Please provide the SeqScreen workflows directory with --workflows_dir (contains .nf workflow files)"
    }
    if (!params.module_dir) {
        error "Please provide the SeqScreen modules directory with --module_dir (contains .sh scripts)"
    }

    // Create input channels from directory containing FASTA and FASTQ files
    // Support multiple file formats: .fa, .fasta, .fas, .fna, .fq, .fastq and their gzipped versions
    
    // Channel for FASTA files (process directly)
    ch_fasta = Channel
        .fromPath("${params.input_dir}/*.{fa,fasta,fas,fna,fa.gz,fasta.gz,fas.gz,fna.gz}")
        .map { fasta ->
            def meta = [
                id: getFileId(fasta),
                single_end: true,
                file_type: 'fasta'
            ]
            tuple(meta, fasta)
        }

    // Channel for FASTQ files (need conversion)
    ch_fastq = Channel
        .fromPath("${params.input_dir}/*.{fq,fastq,fq.gz,fastq.gz}")
        .map { fastq ->
            def meta = [
                id: getFileId(fastq),
                single_end: true,
                file_type: 'fastq'
            ]
            tuple(meta, fastq)
        }

    // Convert FASTQ files to FASTA using seqkit fq2fa
    SEQKIT_FQ2FA(ch_fastq)
    ch_converted_fasta = SEQKIT_FQ2FA.out.fasta

    // Combine FASTA files (original + converted from FASTQ)
    ch_all_fasta = ch_fasta.mix(ch_converted_fasta)

    // Create database channel
    ch_databases = params.databases

    // Create configuration map (added new rapsearch2 parameters)
    ch_config = [
        mode: params.mode,
        threads: params.threads,
        working: params.working,
        evalue: params.evalue,
        slurm: params.slurm,
        version: params.version,
        prefix: params.prefix,
        workflows_dir: params.workflows_dir,
        bin_dir: params.bin_dir,
        module_dir: params.module_dir,
        log: params.log,
        hmmscan: params.hmmscan,
        blastn: params.blastn,
        ancestral: params.ancestral,
        includecc: params.includecc,
        bitscore_cutoff: params.bitscore_cutoff,
        taxlimit: params.taxlimit,
        splitby: params.splitby,
        max_rapsearch_threads: params.max_rapsearch_threads,
        min_file_size_mb: params.min_file_size_mb
    ]

    // Run the SeqScreen workflow
    SEQSCREEN(
        ch_all_fasta,
        ch_databases,
        ch_config
    )

    // Optional: View results
    SEQSCREEN.out.report_html_ch.view { meta, html -> "HTML report for ${meta.id}: ${html}" }
    SEQSCREEN.out.tables_ch.view { meta, table -> "Table report for ${meta.id}: ${table}" }
}