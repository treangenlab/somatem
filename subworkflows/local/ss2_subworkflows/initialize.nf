#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.help = false
params.fasta = ''
params.working = ''
params.log = '/dev/null'
params.slurm = false
params.sensitive = ''
params.bin_dir = null
params.module_dir = null

mode = "fast"
Executor = 'local'

def usage() {
    log.info ''
    log.info 'Usage: nextflow run initialize.nf --fasta /Path/to/infile.fasta --working /Path/to/working_directory [--slurm] [--log=/Path/to/run.log]'
    log.info '  --fasta            Path to the input NT FASTA file'
    log.info '  --working          Path to the output working directory'
    log.info '  --bin_dir          Path to SeqScreen bin directory'
    log.info '  --module_dir       Path to SeqScreen modules directory'
    log.info '  --slurm            Submit modules in this workflow to run on a SLURM grid (Default = run locally)'
    log.info "  --log              Log file that the status can be tee'd to (Default=Don't save the log)"
    log.info '  --help             Print this help message out'
    log.info ''
    exit 1
}

if (params.help) {
    usage()
}

if (params.slurm) {
    Executor = 'slurm'
}

if (!params.fasta) {
    log.error "Missing argument for --fasta option"
    usage()
}

if (!params.working) {
    log.error "Missing argument for --working option"
    usage()
}

logfile = file(params.log)
fastafile = file(params.fasta)
workingDir = file(params.working)

BASE = fastafile.getName()
MODULES = "${params.module_dir}"
SCRIPTS = "${params.bin_dir}"
WORKFLOW = "initialize"

process CREATE_WORKING_DIRECTORIES {
    conda 'bioconda::seqscreen'
    
    output:
    val true, emit: working_dir

    script:
    """
    mkdir -p ${workingDir}
    mkdir -p ${workingDir}/${WORKFLOW}
    """
}

process VERIFY_FASTA {
    conda 'bioconda::seqscreen'
    
    input:
    val working_dir

    output:
    val true, emit: fasta_verified

    executor Executor

    script:
    """
    echo -n " # Launching ${WORKFLOW} workflow .................... " | tee -a ${logfile}; date '+%H:%M:%S %Y-%m-%d' | tee -a ${logfile}

    # Run FASTA validation using seqscreen tools
    ${SCRIPTS}/validate_fasta.pl -f ${fastafile} --max_seq_size=1000000000
    """
}

process INITIALIZE {
    conda 'bioconda::seqscreen'
    
    input:
    val fasta_verified

    output:
    val true, emit: initialized

    executor Executor

    script:
    """
    echo -n " # ${WORKFLOW} workflow complete ......................... " | tee -a ${logfile}; date '+%H:%M:%S %Y-%m-%d' | tee -a ${logfile}
    
    # Create a simple completion marker
    echo "INITIALIZE_COMPLETE=true" > ${workingDir}/initialize_complete.txt
    """
}

workflow {
    CREATE_WORKING_DIRECTORIES()
    VERIFY_FASTA(CREATE_WORKING_DIRECTORIES.out)
    INITIALIZE(VERIFY_FASTA.out)
}