#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.fasta = ''
params.working = ''
params.databases = ''
params.threads = 1
params.help = false
params.sensitive = false
params.ont = false
params.slurm = false
params.hmmscan = false
params.version = false
params.ancestral = false
params.report_only = false
params.skip_report = false
params.check_install = false
params.blastn = false
params.taxlimit = 25
params.report_prefix = false
params.includecc = false
params.run_unmapped = false
params.online = false
params.bitscore = -1
params.splitby = 100000
params.format = 1
params.taxonomy_confidence_threshold = 0
params.evalue = ""
params.filter_taxon = "\"\""
params.keep_taxon = "\"\""
params.keep_html_ont = false
params.window = false
params.window_length = 200
params.window_overlap = 100
slurm_flag = ''
mode = "fast"
sensitive_flag = ''
evalue_flag = ''
hmmscan_flag = ''
ancestral_flag = ''
prefix_flag = ''
include_cc_flag = ""
blastn_flag = ""
ont_flag = ""
window_flag = ""
online_flag = ""
filter_taxon_flag = "--filter-taxon \"\""
keep_taxon_flag = "--keep-taxon \"\""
keep_html_ont_flag = ""
skip_report_flag = ""
VERSION="4.5"

Executor = "local"

WORKFLOWS =  "$workflow.projectDir/workflows"
SCRIPTS = "$workflow.projectDir/scripts"

def usage(ret=1) {
    log.info ""
    log.info "Usage: seqscreen2 --fasta /Path/to/infile.fasta --databases /Path/to/databases --working /Path/to/working_directory [OPTIONS]..."
    log.info "  --fasta         Path to the input NT FASTA file"
    log.info "  --databases     Path to the databases directory containing centrifuge, blast, go etc"
    log.info "  --working       Path to the output working directory"
    log.info "  --threads       Number of threads to use (Default=1)"
    log.info "  --sensitive     Use SeqScreen sensitive mode (old default mode)"
    log.info "  --ont           Enable SeqScreen to run ONT reads"
    log.info "  --window        Run SeqScreen in windowed mode, i.e. run on individual pieces of input sequences of specified length and overlap"
    log.info "  --window_length In windowed mode, use sequence windows of this size. (Default = 200)"
    log.info "  --window_overlap    For windowed mode, overlap between windows (Default = 100)"
    log.info "  --evalue        Cutoff to use for blastx/diamond"
    log.info '  --hmmscan       Run hmmscan on input sequences'
    log.info "  --bitscore      Tiebreak across all proteins within this % of the top bitscore" 
    log.info '  --ancestral     Include all ancestral GO terms in output'
    log.info "  --splitby       Max number of sequences in an input chunk to diamond"
    log.info "  --includecc     Include cellular component go terms"
    log.info "  --blastn        Run blastn in addition to blastx (sensitive mode only)"
    log.info "  --taxlimit      Maximum number of multi-taxIDs to output for a single query in fast mode"
    log.info "  --slurm         Have pipeline modules run on SLURM execution nodes (Default = run locally)"
    log.info "  --report_prefix Add prefix to beginning of seqscreen_report.tsv and seqscreen_html_report.zip. The prefix will be the basename of the input file"
    log.info "  --skip_report   Skip report generation step and only generate intermediate files" 
    log.info "  --report_only   Remove intermediate output and only save the results in {output_dir}/report_generation"
    log.info "  --format        Format type: [1] Original, [2] Hits only, [3] FunSoC only, [4] Gene-Centric,[5] Gene-Centric FunSoC Only [Default: 1]"
    log.info "  --online        Pull reference genomes from NCBI for reference_inference [Needs web access]"
    log.info "  --filter_taxon  Filter comma separated list of taxon"
    log.info "  --keep_taxon    Keep comma separated list of taxon"
    log.info "  --taxonomy_confidence_threshold Confidence threshold for multi-taxids (Average) [Default 0.0]"
    log.info "  --keep_html_ont Keep html report in ont mode [Takes additional memory]" 
    log.info "  --check_install Check for required command line tools, python imports, and database files"
    log.info "  --version       Display the version and exit"
    log.info "  -h, --help      Print this help message out"
    log.info ""
    exit ret
}

if (params.help) {
    usage(0)
}

if (params.version) {
    log.info "SeqScreen v${VERSION}"
    exit 0
}

if (params.bitscore >= 0) {
    bitscore_flag = "--bitscore_cutoff ${params.bitscore}"
} else {
    bitscore_flag = ""
}

if (!params.databases) {
    log.error "Missing argument for --databases option"
    usage()
}

if (params.includecc) {
    include_cc_flag = "--include-cc"
}

if (params.blastn) {
    blastn_flag = "--blastn"
}

if (params.slurm) {
   slurm_flag = '--slurm'
}

if (!params.fasta) {
    log.error "Missing argument for --fasta option"
    usage()
}

if (!params.working) {
    log.error "Missing argument for --working option"
    usage()
}

if (params.sensitive) {
    mode = "sensitive"
    sensitive_flag = "--sensitive"
}

if (params.ont) {
    mode = "ont"
}

if (params.hmmscan) {
    hmmscan_flag = "--hmmscan"
}
if (params.ancestral) {
    ancestral_flag = "--ancestral"
}

if (params.skip_report) {
    skip_report_flag = "--skip_report"
}

if (params.evalue != "") {
    evalue_flag = "--evalue " + params.evalue
}
if (params.online) {
   online_flag = "--online"
}

if (params.filter_taxon) {
    if (mode == "ont") {
        filter_taxon_flag = "--filter_taxon ${params.filter_taxon}"
    }
    else {
        filter_taxon_flag = "--filter-taxon ${params.filter_taxon}"
    }
}

if (params.keep_taxon) {
    if (mode == "ont") {
        keep_taxon_flag = "--keep_taxon ${params.keep_taxon}"
    }
    else {
        keep_taxon_flag = "--keep-taxon ${params.keep_taxon}"
    }
} 

if (params.keep_html_ont) {
   keep_html_ont_flag = "--keep_html_ont"
}

fastafile = file(params.fasta)
workingDir = file(params.working)
databasesDir = file(params.databases)
taxDir = file("${workingDir}/taxonomic_identification/taxonomic_assignment")
logfile = "$workingDir/seqscreen.log"
BASE = fastafile.getName()

if (params.check_install) {
    result = "python ${SCRIPTS}/check_install.py -d ${databasesDir}".execute().text
    print(result)
    exit 0
}

if (params.window) {
    window_flag = "--window"
}
//     newfastafile = "${workingDir}/${BASE}_windowed.fasta"
// } else {
//     newfastafile = fastafile
// }

if (params.report_prefix) {
    prefix_flag = "--prefix ${BASE}_"
}

process VERIFY_FASTA {
    
    publishDir "${workingDir}", mode: 'copy', pattern: 'seqscreen.log'

    output:
        path "seqscreen.log", emit: log_file

    script:
    """
    mkdir -p ${workingDir}
    if [ -e seqscreen.log ]; then rm seqscreen.log; fi    
    echo -n " #### Launching SeqScreen pipeline version ${VERSION} ....... " | tee -a seqscreen.log
    date '+%H:%M:%S %Y-%m-%d' | tee -a seqscreen.log
    """
}

process WINDOW_FASTA {

    input:
        path fasta

    output:
        path "*.fasta", emit: windowed

    publishDir "${workingDir}", mode: 'copy'

    script:
    """
    base=\$(basename ${fasta})
    base=\${base%.*}
    output_file="\${base}_windowed.fasta"

    python3 ${SCRIPTS}/split_window_input.py --fasta ${fasta} \
                                             --output "\${output_file}" \
                                             --length ${params.window_length} \
                                             --overlap ${params.window_overlap}
    """
}

process INITIALIZE {

    input:
        path fasta

    output:
        val true, emit: initialized

    script:
    """
    nextflow run ${WORKFLOWS}/initialize.nf --fasta ${fasta}\
                                                 --databases ${databasesDir} \
                                                 --working ${workingDir} \
                                                 ${params.slurm ? '--slurm' : ''} \
                                                 ${params.sensitive ? '--sensitive' : ''} \
                                                 --log ${logfile}
    """
}

process ONT {

    input:
        val initialized
        path fasta

    output:
        val true, emit: ont

    script:
    """
    nextflow run ${WORKFLOWS}/nanopore.nf --fasta ${fasta} \
                                        --databases ${databasesDir} \
                                        --working ${workingDir} \
                                        --threads ${params.threads} \
                                        --splitby ${params.splitby} \
                                        $prefix_flag \
                                        $bitscore_flag \
                                        $ancestral_flag \
                                        $include_cc_flag \
                                        $evalue_flag \
                                        $slurm_flag \
                                        $skip_report_flag \
                                        --format ${params.format} \
                                        $online_flag \
                                        $filter_taxon_flag \
                                        $keep_taxon_flag \
                                        $keep_html_ont_flag \
                                        --taxonomy_confidence_threshold ${params.taxonomy_confidence_threshold} \
                                        --version $VERSION \
                                        --log ${logfile}
    """
}

process SEQMAPPER {

    input:
        val initialized
        path fasta

    output:
        val true, emit: seqmapper

    script:
    """
    nextflow run ${WORKFLOWS}/seqmapper.nf --fasta ${fasta} \
                                            --databases ${databasesDir} \
                                            --working ${workingDir} \
                                            --threads ${params.threads} \
                                            $slurm_flag \
                                            $hmmscan_flag \
                                            --log ${logfile}
    """
}

process TAXONOMIC_IDENTIFICATION {

    input:
        val initialized
        path fasta

    output:
        val true, emit: taxonomic_identification

    script:
    if (mode == "fast") {
        """
        nextflow run ${WORKFLOWS}/taxonomic_identification_fast.nf --fasta ${fasta} \
                                            --databases ${databasesDir} \
                                            --working ${workingDir} \
                                            --threads ${params.threads} \
                                                                   --splitby ${params.splitby} \
                                                                   --taxlimit ${params.taxlimit} \
                                                                   $slurm_flag \
                                                                   $evalue_flag \
                                                                   --log ${logfile}                                                         
        """
    }
    else if (mode == "sensitive") {
        """
        nextflow run ${WORKFLOWS}/taxonomic_identification_sensitive.nf --fasta ${fasta} \
                                            --databases ${databasesDir} \
                                            --working ${workingDir} \
                                            --threads ${params.threads} \
                                                                       $slurm_flag \
                                                                       $blastn_flag \
                                                                       $evalue_flag \
                                                                       --log ${logfile}                                                             
        """
    }
}

process FUNCTIONAL_ANNOTATION {

    input:
        val taxonomic_identification
        path fasta

    output:
        val true, emit: functional_annotation

    script:
    """
    nextflow run ${WORKFLOWS}/functional_annotation.nf --fasta ${fasta} \
                                                     --databases ${databasesDir} \
                                                     --working ${workingDir} \
                                                     --threads ${params.threads} \
                                                     $bitscore_flag \
                                                     $ancestral_flag \
                                                     $include_cc_flag \
                                                     $slurm_flag \
                                                     $sensitive_flag \
                                                     --log ${logfile}                                            
    """
}

process REPORT_GENERATION {

    input:
        val functional_annotation
        val taxonomic_identification
        path fasta

    output:
        val true, emit: report_generation

    script:
    """
    nextflow run ${WORKFLOWS}/report_generation.nf --working ${workingDir} \
                                                        --version ${VERSION} \
                                                        --fasta ${fasta} \
                                                        --databases ${databasesDir} \
                                                        $prefix_flag \
                                                        $slurm_flag \
                                                        $ancestral_flag \
                                                        $hmmscan_flag \
                                                        $evalue_flag \
                                                        --bitscore ${params.bitscore} \
                                                        $blastn_flag \
                                                        --splitby ${params.splitby} \
                                                        --taxlimit ${params.taxlimit} \
                                                        --includecc ${params.includecc} \
                                                        $sensitive_flag \
                                                        --log ${logfile}
    echo -n " #### SeqScreen pipeline complete ......................... " | tee -a ${logfile}; date '+%H:%M:%S %Y-%m-%d' | tee -a ${logfile}
    """
}

process INFERENCE {

    input:
        val report_generation
        path fasta

    output:
        val true, emit: inference
    
    executor Executor

    script:
    """
    echo -n "# Launching inference workflow ...... " | tee -a ${logfile}; date '+%H:%M:%S %Y-%m-%d' | tee -a ${logfile}
    python3 ${SCRIPTS}/reference_inference_short_read.py --fasta1 ${fasta} \
                                              --output ${workingDir} \
                                              --working ${workingDir}/taxonomic_identification/taxonomic_assignment/inference_working \
                                              --databases ${databasesDir} \
                                              --threads ${params.threads} \
                                              $online_flag
    echo -n "# inference complete ......." | tee -a ${logfile}; date '+%H:%M:%S %Y-%m-%d' | tee -a ${logfile}
    """
}

process FORMAT {

    input:
        val inference

    output:
        val true, emit: formatted

    executor Executor

    script:
    if (params.skip_report)
   """
   """
   else
   """
   echo -n "# Launching format workflow ...... " | tee -a ${logfile}; date '+%H:%M:%S %Y-%m-%d' | tee -a ${logfile}
   python3 ${SCRIPTS}/format.py --report ${workingDir}/report_generation \
                                --format ${params.format} \
                                --mode ${mode} \
                                --databases ${databasesDir} \
                                --taxonomy_confidence_threshold ${params.taxonomy_confidence_threshold} \
                                $filter_taxon_flag \
                                $keep_taxon_flag 
   """
}

process CLEANUP {

    input:
        val formatted

    output:
        val true, emit: cleanup

    script:
    if (params.sensitive) {
        """
        echo -n "# Launching cleanup workflow ...... " | tee -a ${logfile}; date '+%H:%M:%S %Y-%m-%d' | tee -a ${logfile}
        rm -r ${workingDir}/functional_annotation ${workingDir}/taxonomic_identification ${workingDir}/seqmapper
        """
    } else {
        """
        echo -n "# Launching cleanup workflow ...... " | tee -a ${logfile}; date '+%H:%M:%S %Y-%m-%d' | tee -a ${logfile}
        rm -r ${workingDir}/functional_annotation ${workingDir}/taxonomic_identification 
        """
    }
}


workflow {
    fasta_ch = Channel.fromPath(params.fasta)

    VERIFY_FASTA()

    if(params.window) {
        WINDOW_FASTA(fasta_ch)
        fasta_for_processing = WINDOW_FASTA.out.windowed
    } else {
        fasta_for_processing = fasta_ch
    }

    INITIALIZE(fasta_for_processing)

    if (params.ont) {
        ONT(INITIALIZE.out, fasta_for_processing)
    } else {
        SEQMAPPER(INITIALIZE.out, fasta_for_processing)
        TAXONOMIC_IDENTIFICATION(INITIALIZE.out, fasta_for_processing)
        FUNCTIONAL_ANNOTATION(TAXONOMIC_IDENTIFICATION.out, fasta_for_processing)
        REPORT_GENERATION(FUNCTIONAL_ANNOTATION.out, TAXONOMIC_IDENTIFICATION.out, fasta_for_processing)
        INFERENCE(REPORT_GENERATION.out, fasta_for_processing)
        FORMAT(INFERENCE.out)
        if (params.report_only) {
            CLEANUP(FORMAT.out)
        }
    }
}