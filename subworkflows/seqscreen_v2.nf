#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
SeqScreen that composes the repo's .nf workflows:

Initialization
SeqMapper
Nanopore (optional - if --sequencing_type long_read)
Taxonomic Identification (Fast or Sensitive mode based on cfg.mode parameter)
Functional Annotation
Report Generation

Usage: nextflow run seqscreen.nf
    --input_dir /path/to/input/directory
    --databases /path/to/seqscreen/databases
    --bin_dir /path/to/seqscreen/bin
    --workflows_dir /path/to/seqscreen/workflows
    --module_dir /path/to/seqscreen/modules
*/

// Default parameter values to avoid undefined parameter warnings
params.mode = 'fast'
params.sequencing_type = 'short_read'
params.threads = 1
params.working = 'work'
params.evalue = 10
params.slurm = false
params.version = '1.0'
params.prefix = ''
params.max_rapsearch_threads = 16
params.min_file_size_mb = 1
params.splitby = 100000
params.taxonomy_confidence_threshold = 0
params.bitscore_cutoff = 5
params.ancestral = false
params.includecc = false
params.skip_report = false
params.online = false
params.format = 1
params.keep_html_ont = false
params.filter_taxon = ''
params.keep_taxon = ''
params.log = '/dev/null'
params.hmmscan = false
params.blastn = false
params.taxlimit = 25
params.help = false

// Help message
def helpMessage() {
    log.info"""
    Usage: nextflow run seqscreen.nf [options]

    Required Arguments:
      --input_dir             Path to input directory containing FASTA files
      --databases             Path to SeqScreen database directory
      --bin_dir               Path to SeqScreen bin directory
      --workflows_dir         Path to SeqScreen workflows directory
      --module_dir            Path to SeqScreen modules directory

    Optional Arguments:
      --working               Working directory (default: work)
      --mode                  Analysis mode: fast or sensitive (default: fast)
      --sequencing_type       Sequencing type: short_read or long_read (default: short_read)
      --threads               Number of threads (default: 1)
      --evalue                E-value threshold (default: 10)
      --slurm                 Enable SLURM execution (default: false)
      --version               SeqScreen version (default: 1.0)
      --prefix                Output file prefix (default: "")
      --max_rapsearch_threads Maximum threads for rapsearch2 (default: 16)
      --min_file_size_mb      Minimum file size in MB for full resources (default: 1)
      --splitby               Number of fasta sequences in each chunk for nanopore (default: 100000)
      --taxonomy_confidence_threshold Confidence threshold for multi-tax ids (default: 0)
      --bitscore_cutoff       Bitscore cutoff for nanopore analysis (default: 5)
      --ancestral             Include all ancestral GO terms (nanopore)
      --includecc             Include cellular component GO terms (nanopore)
      --skip_report           Skip report generation step (nanopore)
      --online                Pull reference genomes from NCBI (nanopore)
      --format                Format type for nanopore output (default: 1)
      --keep_html_ont         Keep html report in ont mode (nanopore)
      --filter_taxon          Filter comma separated list of taxon (nanopore)
      --keep_taxon            Keep comma separated list of taxon (nanopore)

    Supported File Formats:
      FASTA: .fa, .fasta, .fas, .fna (and their .gz versions)

    Profiles:
      -profile standard       Run locally (default)
      -profile slurm          Run on SLURM cluster

    Example:
      nextflow run seqscreen.nf \\
        --input_dir /path/to/sequences/ \\
        --databases /path/to/db \\
        --workflows_dir /path/to/workflows \\
        --bin_dir /path/to/bin \\
        --module_dir /path/to/modules \\
        --sequencing_type long_read \\
        --threads 8
    """.stripIndent()
}

// Show help message if requested
if (params.help) {
    helpMessage()
    exit 0
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
    val db_path_str
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

    # Debug database path
    echo "=== DATABASE PATH DEBUG ===" | tee -a ${log_file}
    echo "Database path received: ${db_path_str}" | tee -a ${log_file}
    echo "Database path exists: \$(test -d '${db_path_str}' && echo 'YES' || echo 'NO')" | tee -a ${log_file}
    if [ -d "${db_path_str}" ]; then
        echo "Database contents:" | tee -a ${log_file}
        ls -la "${db_path_str}/" | head -5 | tee -a ${log_file}
    fi
    echo "=== END DATABASE PATH DEBUG ===" | tee -a ${log_file}

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
    val db_path_str
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
    echo "SeqMapper database path: ${db_path_str}" | tee -a ${log_file}

    nextflow run ${workflows_dir}/seqmapper.nf \\
        --fasta ${fasta} \\
        --working ${working_dir} \\
        --databases "${db_path_str}" \\
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

process RUN_NANOPORE {
    tag "${meta.id}"
    label 'process_high'
    conda 'bioconda::seqscreen'

    input:
    tuple val(meta), path(fasta)
    val db_path_str
    val cfg

    output:
    tuple val(meta), path(fasta), emit: nanopore_results

    script:
    def working_dir = cfg?.working ?: "work"
    def log_file = cfg?.log ?: "/dev/null"
    def threads = cfg?.threads ?: 1
    def evalue = cfg?.evalue ?: 0.001
    def splitby = cfg?.splitby ?: 100000
    def taxonomy_confidence_threshold = cfg?.taxonomy_confidence_threshold ?: 0
    def bitscore_cutoff = cfg?.bitscore_cutoff ?: 5
    def slurm_flag = cfg?.slurm ? "--slurm" : ""
    def ancestral_flag = cfg?.ancestral ? "--ancestral" : ""
    def includecc_flag = cfg?.includecc ? "--includecc" : ""
    def skip_report_flag = cfg?.skip_report ? "--skip_report" : ""
    def online_flag = cfg?.online ? "--online" : ""
    def keep_html_ont_flag = cfg?.keep_html_ont ? "--keep_html_ont" : ""
    def filter_taxon = cfg?.filter_taxon ?: ""
    def keep_taxon = cfg?.keep_taxon ?: ""
    def format = cfg?.format ?: 1
    def prefix = cfg?.prefix ?: ""
    def version = cfg?.version ?: "1.0"
    def workflows_dir = cfg?.workflows_dir ?: "."
    def bin_dir = cfg?.bin_dir ?: "."
    def module_dir = cfg?.module_dir ?: "."

    """
    # Set up parameters and run nanopore workflow
    export NXF_WORK=\${PWD}/work_nanopore

    # Set environment variables that the sub-subworkflow expects
    export BIN_DIR=${bin_dir}
    export MODULE_DIR=${module_dir}

    echo "Nanopore database path: ${db_path_str}" | tee -a ${log_file}

    nextflow run ${workflows_dir}/nanopore.nf \\
        --fasta ${fasta} \\
        --working ${working_dir} \\
        --databases "${db_path_str}" \\
        --bin_dir ${bin_dir} \\
        --module_dir ${module_dir} \\
        --threads ${threads} \\
        --evalue ${evalue} \\
        --splitby ${splitby} \\
        --taxonomy_confidence_threshold ${taxonomy_confidence_threshold} \\
        --bitscore_cutoff ${bitscore_cutoff} \\
        --format ${format} \\
        --prefix ${prefix} \\
        --version ${version} \\
        --filter_taxon "${filter_taxon}" \\
        --keep_taxon "${keep_taxon}" \\
        ${slurm_flag} \\
        ${ancestral_flag} \\
        ${includecc_flag} \\
        ${skip_report_flag} \\
        ${online_flag} \\
        ${keep_html_ont_flag} \\
        -work-dir \${NXF_WORK}
    """
}

process RUN_TAXONOMIC_IDENTIFICATION_FAST {
    tag "${meta.id}"
    label 'process_high'
    conda 'bioconda::seqscreen'

    input:
    tuple val(meta), path(fasta)
    val db_path_str
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

    echo "Taxonomic ID Fast database path: ${db_path_str}" | tee -a ${log_file}

    nextflow run ${workflows_dir}/taxonomic_identification_fast.nf \\
        --fasta ${fasta} \\
        --working ${working_dir} \\
        --databases "${db_path_str}" \\
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
    val db_path_str
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

    echo "Taxonomic ID Sensitive database path: ${db_path_str}" | tee -a ${log_file}

    nextflow run ${workflows_dir}/taxonomic_identification_sensitive.nf \\
        --fasta ${fasta} \\
        --working ${working_dir} \\
        --databases "${db_path_str}" \\
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
    val db_path_str
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

    echo "Functional annotation database path: ${db_path_str}" | tee -a ${log_file}

    nextflow run ${workflows_dir}/functional_annotation.nf \\
        --fasta ${fasta} \\
        --working ${working_dir} \\
        --databases "${db_path_str}" \\
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
    val nanopore_complete         // Signal that nanopore is done (if applicable)
    val db_path_str
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

    echo "Report generation database path: ${db_path_str}" | tee -a ${log_file}

    nextflow run ${workflows_dir}/report_generation.nf \\
        --fasta ${fasta} \\
        --working ${working_dir} \\
        --databases "${db_path_str}" \\
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
        echo '{"status": "completed", "sample": "'${meta.id}'", "sequencing_type": "'${cfg?.sequencing_type ?: 'short_read'}'"}' > report_output/summary.json
    fi
    """
}

workflow SEQSCREEN {

    take:
    seqs_ch
    db_path_str
    cfg

    main:
    //
    // Adapter layer
    // Make tiny, explicit tuples so the ports line up regardless of your upstream.
    //
    ch_input = seqs_ch
    ch_db_str = db_path_str
    ch_cfg   = cfg

    //
    // 1) Initialization
    // Initialize the SeqScreen workflow and validate input
    //
    RUN_INITIALIZE(ch_input, ch_db_str, ch_cfg)
    init_seqs = RUN_INITIALIZE.out.initialized

    //
    // 2) Mapping (ALWAYS RUN)
    // Run sequence mapping against databases
    //
    RUN_SEQMAPPER(init_seqs, ch_db_str, ch_cfg)
    mapped = RUN_SEQMAPPER.out.mapped

    //
    // 3) Nanopore Analysis (CONDITIONAL)
    // Run nanopore workflow if sequencing_type is long_read
    //
    if (ch_cfg.sequencing_type == "long_read") {
        RUN_NANOPORE(mapped, ch_db_str, ch_cfg)
        nanopore_results = RUN_NANOPORE.out.nanopore_results
        nanopore_complete = nanopore_results.map { it -> true }
    } else {
        // Create a dummy completion signal for short reads
        nanopore_complete = Channel.value(true)
    }

    //
    // 4) Taxonomic identification
    // Call different sub-subworkflows based on user parameters (mode: fast/sensitive)
    // IMPORTANT: This must complete before functional annotation starts
    //

    // Determine which taxonomic identification workflow to run based on mode
    if (ch_cfg.mode == "sensitive") {
        // Run sensitive taxonomic identification
        RUN_TAXONOMIC_IDENTIFICATION_SENSITIVE(
            mapped,
            ch_db_str,
            ch_cfg
        )
        tax_calls = RUN_TAXONOMIC_IDENTIFICATION_SENSITIVE.out.taxonomic_results
    } else {
        // Run fast taxonomic identification (default)
        RUN_TAXONOMIC_IDENTIFICATION_FAST(
            mapped,
            ch_db_str,
            ch_cfg
        )
        tax_calls = RUN_TAXONOMIC_IDENTIFICATION_FAST.out.taxonomic_results
    }

    //
    // 5) Functional annotation
    // FIXED: Now takes input from taxonomic identification (tax_calls) instead of mapped
    // This ensures taxonomic identification completes before functional annotation starts
    //
    RUN_FUNCTIONAL_ANNOTATION(tax_calls, ch_db_str, ch_cfg)
    func_calls = RUN_FUNCTIONAL_ANNOTATION.out.functional_results

    //
    // 6) Report generation
    // FIXED: Now takes only tax_calls (with FASTA) and completion signals from functional annotation and nanopore
    // This avoids the file name collision while ensuring both processes complete first
    //
    RUN_REPORT_GENERATION(
        tax_calls,                    // Contains meta and fasta
        func_calls.map { it -> true }, // Just a completion signal, no files
        nanopore_complete,            // Completion signal from nanopore (or dummy for short reads)
        ch_db_str, 
        ch_cfg
    )

    //
    // 7) Emit canonical ports you can hook into MultiQC/Dash/etc.
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

    // Debug database path at the very beginning
    log.info "=== MAIN WORKFLOW DATABASE DEBUG ==="
    log.info "Database parameter received: ${params.databases}"
    log.info "Database path type: ${params.databases.getClass()}"
    
    // Create input channels from directory containing FASTA files only
    // Support multiple FASTA file formats: .fa, .fasta, .fas, .fna and their gzipped versions
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

    // Create database channel as a simple string to avoid path corruption
    ch_databases = params.databases.toString()

    // Create configuration map (added nanopore-specific parameters)
    ch_config = [
        mode: params.mode,
        sequencing_type: params.sequencing_type,
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
        min_file_size_mb: params.min_file_size_mb,
        // Nanopore-specific parameters
        taxonomy_confidence_threshold: params.taxonomy_confidence_threshold,
        skip_report: params.skip_report,
        online: params.online,
        format: params.format,
        keep_html_ont: params.keep_html_ont,
        filter_taxon: params.filter_taxon,
        keep_taxon: params.keep_taxon
    ]

    // Run the SeqScreen workflow
    SEQSCREEN(
        ch_fasta,
        ch_databases,
        ch_config
    )

    // Optional: View results
    SEQSCREEN.out.report_html_ch.view { meta, html -> "HTML report for ${meta.id}: ${html}" }
    SEQSCREEN.out.tables_ch.view { meta, table -> "Table report for ${meta.id}: ${table}" }
}