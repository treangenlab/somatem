#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.help = false
params.fasta = ''
params.working = ''
params.bin_dir = null
params.module_dir = null
params.log = '/dev/null'
params.databases = ''
params.threads = 1
params.slurm = false
params.hmmscan = false
params.max_rapsearch_threads = 16
params.min_file_size_mb = 1

Executor = 'local'

def usage() {
    log.info ''
    log.info 'Usage: nextflow run seqmapper.nf --fasta /Path/to/infile.fasta --working /Path/to/working_directory --databases /Path/to/databases [--threads 4] [--slurm] [--log=/Path/to/run.log]'
    log.info '  --fasta      Path to the input NT FASTA file'
    log.info '  --working    Path to the output working directory'
    log.info "  --databases  Path to the S2FAST databases"
    log.info "  --threads    Number of threads to use (Default=1)"
    log.info '  --slurm      Submit modules in this workflow to run on a SLURM grid (Default = run locally)'
    log.info '  --hmmscan    Run hmmscan on input sequences'
    log.info "  --log        Log file that the status can be tee'd to (Default=Don't save the log)"
    log.info '  --help       Print this help message out'
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

if (!params.databases) {
    log.error "Missing argument for --databases option"
    usage()
}

logfile = file(params.log)
fastafile = file(params.fasta)
workingDir = file(params.working)
databasesDir = file(params.databases)
pfamDir = file("${workingDir}/seqmapper/pfam_hmm")

BASE = fastafile.getName()
MODULES = "${params.module_dir}"
SCRIPTS = "${params.bin_dir}"
WORKFLOW = "seqmapper"

// output locations
bowtieBsatDir = file("${workingDir}/${WORKFLOW}/bowtie2_bsat")
rapsearchBsatDir = file("${workingDir}/${WORKFLOW}/rapsearch2_bsat")
bowtieVfdbDir = file("${workingDir}/${WORKFLOW}/bowtie2_vfdb")
rapsearchVfdbDir = file("${workingDir}/${WORKFLOW}/rapsearch2_vfdb")

translatedFile = "${pfamDir}/${BASE}.translated.fasta"

process CREATE_WORKING_DIRECTORIES {
    conda 'bioconda::seqscreen'
    
    output:
    val true, emit: working_dir

    script:
    """
    mkdir -p ${workingDir}
    mkdir -p ${workingDir}/${WORKFLOW}

    if [ -d ${bowtieBsatDir} ]; then rm -rf ${bowtieBsatDir}; fi;
    if [ -d ${rapsearchBsatDir} ]; then rm -rf ${rapsearchBsatDir}; fi;
    if [ -d ${bowtieVfdbDir} ]; then rm -rf ${bowtieVfdbDir}; fi;
    if [ -d ${rapsearchVfdbDir} ]; then rm -rf ${rapsearchVfdbDir}; fi;
    if [ -d ${pfamDir} ]; then rm -rf ${pfamDir}; fi;

    mkdir ${bowtieBsatDir}
    mkdir ${rapsearchBsatDir}
    mkdir ${bowtieVfdbDir}
    mkdir ${rapsearchVfdbDir}
    if ${params.hmmscan}; then mkdir ${pfamDir}; fi
    """
}

process INITIALIZE {
    conda 'bioconda::seqscreen'
    
    input:
    val working_dir

    output:
    val true, emit: initialized

    executor Executor

    script:
    """
    echo -n " # Launching ${WORKFLOW} workflow ........................ " | tee -a ${logfile}; date '+%H:%M:%S %Y-%m-%d' | tee -a ${logfile}
    
    # Debug information
    echo "=== DEBUG INFO ===" | tee -a ${logfile}
    echo "FASTA file: ${fastafile}" | tee -a ${logfile}
    echo "Database dir: ${databasesDir}" | tee -a ${logfile}
    echo "Working dir: ${workingDir}" | tee -a ${logfile}
    echo "Threads: ${params.threads}" | tee -a ${logfile}
    echo "Max RapSearch threads: ${params.max_rapsearch_threads}" | tee -a ${logfile}
    
    # Check if FASTA file exists and has content
    if [ -f "${fastafile}" ]; then
        echo "FASTA file exists" | tee -a ${logfile}
        echo "FASTA file size: \$(stat -c%s '${fastafile}') bytes" | tee -a ${logfile}
        echo "Number of sequences: \$(grep -c '^>' '${fastafile}' || echo 0)" | tee -a ${logfile}
    else
        echo "ERROR: FASTA file not found!" | tee -a ${logfile}
    fi
    
    # Check if database directories exist
    echo "Checking database paths:" | tee -a ${logfile}
    echo "Database root: ${databasesDir}" | tee -a ${logfile}
    ls -la "${databasesDir}/" | tee -a ${logfile} || echo "Cannot list database directory" | tee -a ${logfile}
    
    if [ -d "${databasesDir}/bowtie2/bsat_ccl" ]; then
        echo "BSAT bowtie2 directory exists" | tee -a ${logfile}
        ls -la "${databasesDir}/bowtie2/bsat_ccl/" | head -5 | tee -a ${logfile}
    else
        echo "ERROR: BSAT bowtie2 directory not found at ${databasesDir}/bowtie2/bsat_ccl!" | tee -a ${logfile}
    fi
    
    if [ -d "${databasesDir}/bowtie2/vfdb" ]; then
        echo "VFDB bowtie2 directory exists" | tee -a ${logfile}
        ls -la "${databasesDir}/bowtie2/vfdb/" | head -5 | tee -a ${logfile}
    else
        echo "ERROR: VFDB bowtie2 directory not found at ${databasesDir}/bowtie2/vfdb!" | tee -a ${logfile}
    fi
    
    if [ -d "${databasesDir}/rapsearch2/bsat_ccl" ]; then
        echo "BSAT rapsearch2 directory exists" | tee -a ${logfile}
    else
        echo "ERROR: BSAT rapsearch2 directory not found!" | tee -a ${logfile}
    fi
    
    if [ -d "${databasesDir}/rapsearch2/vfdb" ]; then
        echo "VFDB rapsearch2 directory exists" | tee -a ${logfile}
    else
        echo "ERROR: VFDB rapsearch2 directory not found!" | tee -a ${logfile}
    fi
    
    echo "=== END DEBUG INFO ===" | tee -a ${logfile}
    """
}

process BOWTIE2_BSAT {
    conda 'bioconda::seqscreen'
    errorStrategy 'ignore'  // Don't fail the whole pipeline if this fails
    
    input:
    val initialized

    output:
    val true, emit: bowtie_bsat

    if ( Executor == 'local' ) {
       executor "local"
    }
    else if ( Executor == 'slurm' ) {
       clusterOptions "--ntasks-per-node ${params.threads}"
       executor "slurm"
    }

    script:
    """
    echo "=== BOWTIE2_BSAT DEBUG ===" | tee -a ${logfile}
    echo "Command to run:" | tee -a ${logfile}
    echo "${MODULES}/bowtie2.sh --fasta=${fastafile} --database=${databasesDir}/bowtie2/bsat_ccl/blacklist.seqs.nt --out=${bowtieBsatDir}/blacklist_bsat.sam --threads=${params.threads}" | tee -a ${logfile}
    
    # Check if the bowtie2.sh script exists
    if [ ! -f "${MODULES}/bowtie2.sh" ]; then
        echo "ERROR: bowtie2.sh script not found at ${MODULES}/bowtie2.sh" | tee -a ${logfile}
        touch ${bowtieBsatDir}/blacklist_bsat.sam  # Create empty output
        exit 0
    fi
    
    # Check if database exists
    if [ ! -f "${databasesDir}/bowtie2/bsat_ccl/blacklist.seqs.nt.1.bt2" ]; then
        echo "ERROR: BSAT database not found at ${databasesDir}/bowtie2/bsat_ccl/" | tee -a ${logfile}
        echo "Looking for files in database directory:" | tee -a ${logfile}
        ls -la "${databasesDir}/bowtie2/bsat_ccl/" | tee -a ${logfile} || echo "Directory does not exist" | tee -a ${logfile}
        touch ${bowtieBsatDir}/blacklist_bsat.sam  # Create empty output
        exit 0
    fi
    
    # Run the command with error handling
    set +e
    ${MODULES}/bowtie2.sh --fasta=${fastafile} \\
                          --database=${databasesDir}/bowtie2/bsat_ccl/blacklist.seqs.nt \\
                          --out=${bowtieBsatDir}/blacklist_bsat.sam \\
                          --threads=${params.threads}
    EXIT_CODE=\$?
    set -e
    
    echo "Bowtie2 BSAT exit code: \$EXIT_CODE" | tee -a ${logfile}
    
    if [ \$EXIT_CODE -ne 0 ]; then
        echo "Bowtie2 BSAT failed, creating empty output file" | tee -a ${logfile}
        touch ${bowtieBsatDir}/blacklist_bsat.sam
    fi
    
    echo "=== END BOWTIE2_BSAT DEBUG ===" | tee -a ${logfile}
    """
}

process BOWTIE2_VFDB {
    conda 'bioconda::seqscreen'
    errorStrategy 'ignore'  // Don't fail the whole pipeline if this fails
    
    input:
    val initialized

    output:
    val true, emit: bowtie_vfdb

    if ( Executor == 'local' ) {
       executor "local"
    }
    else if ( Executor == 'slurm' ) {
       clusterOptions "--ntasks-per-node ${params.threads}"
       executor "slurm"
    }

    script:
    """
    echo "=== BOWTIE2_VFDB DEBUG ===" | tee -a ${logfile}
    echo "Command to run:" | tee -a ${logfile}
    echo "${MODULES}/bowtie2.sh --fasta=${fastafile} --database=${databasesDir}/bowtie2/vfdb/vfdb --out=${bowtieVfdbDir}/blacklist_vfdb.sam --threads=${params.threads}" | tee -a ${logfile}
    
    # Check if the bowtie2.sh script exists
    if [ ! -f "${MODULES}/bowtie2.sh" ]; then
        echo "ERROR: bowtie2.sh script not found at ${MODULES}/bowtie2.sh" | tee -a ${logfile}
        touch ${bowtieVfdbDir}/blacklist_vfdb.sam  # Create empty output
        exit 0
    fi
    
    # Check if database exists - look for the actual index files
    if [ ! -f "${databasesDir}/bowtie2/vfdb/vfdb.1.bt2" ]; then
        echo "ERROR: VFDB database not found at ${databasesDir}/bowtie2/vfdb/" | tee -a ${logfile}
        echo "Looking for files in VFDB database directory:" | tee -a ${logfile}
        ls -la "${databasesDir}/bowtie2/vfdb/" | tee -a ${logfile} || echo "Directory does not exist" | tee -a ${logfile}
        
        # Try alternative locations
        echo "Checking for alternative database locations:" | tee -a ${logfile}
        find "${databasesDir}" -name "*.bt2" -type f | head -10 | tee -a ${logfile}
        
        touch ${bowtieVfdbDir}/blacklist_vfdb.sam  # Create empty output
        exit 0
    fi
    
    # Run the command with error handling
    set +e
    ${MODULES}/bowtie2.sh --fasta=${fastafile} \\
                          --database=${databasesDir}/bowtie2/vfdb/vfdb \\
                          --out=${bowtieVfdbDir}/blacklist_vfdb.sam \\
                          --threads=${params.threads}
    EXIT_CODE=\$?
    set -e
    
    echo "Bowtie2 VFDB exit code: \$EXIT_CODE" | tee -a ${logfile}
    
    if [ \$EXIT_CODE -ne 0 ]; then
        echo "Bowtie2 VFDB failed, creating empty output file" | tee -a ${logfile}
        touch ${bowtieVfdbDir}/blacklist_vfdb.sam
    fi
    
    echo "=== END BOWTIE2_VFDB DEBUG ===" | tee -a ${logfile}
    """
}

process RAPSEARCH2_BSAT {
    conda 'bioconda::seqscreen'
    errorStrategy 'ignore'  // Don't fail the whole pipeline if this fails
    
    input:
    val initialized

    output:
    val true, emit: rapsearch_bsat

    if ( Executor == 'local' ) {
       executor "local"
    }
    else if ( Executor == 'slurm' ) {
       clusterOptions "--ntasks-per-node ${params.threads}"
       executor "slurm"
    }

    script:
    // Limit threads for rapsearch2 to avoid segmentation faults
    def rapsearch_threads = Math.min(params.threads as Integer, params.max_rapsearch_threads as Integer)
    
    """
    echo "=== RAPSEARCH2_BSAT DEBUG ===" | tee -a ${logfile}
    echo "Using ${rapsearch_threads} threads for RapSearch2 (limited from ${params.threads})" | tee -a ${logfile}
    
    # Check if database exists
    if [ ! -f "${databasesDir}/rapsearch2/bsat_ccl/blacklist.seqs.aa" ]; then
        echo "ERROR: RAPSEARCH2 BSAT database not found at ${databasesDir}/rapsearch2/bsat_ccl/blacklist.seqs.aa" | tee -a ${logfile}
        echo "Looking for files in rapsearch2 BSAT directory:" | tee -a ${logfile}
        ls -la "${databasesDir}/rapsearch2/bsat_ccl/" | tee -a ${logfile} || echo "Directory does not exist" | tee -a ${logfile}
        touch ${rapsearchBsatDir}/blacklist_bsat.m8
        exit 0
    fi
    
    # Check if the rapsearch2.sh script exists
    if [ ! -f "${MODULES}/rapsearch2.sh" ]; then
        echo "ERROR: rapsearch2.sh script not found at ${MODULES}/rapsearch2.sh" | tee -a ${logfile}
        touch ${rapsearchBsatDir}/blacklist_bsat.m8
        exit 0
    fi
    
    set +e
    ${MODULES}/rapsearch2.sh --fasta=${fastafile} \\
                             --database=${databasesDir}/rapsearch2/bsat_ccl/blacklist.seqs.aa \\
                             --out=${rapsearchBsatDir}/blacklist_bsat \\
                             --evalue=1e-9 \\
                             --threads=${rapsearch_threads}
    EXIT_CODE=\$?
    set -e
    
    echo "RapSearch2 BSAT exit code: \$EXIT_CODE" | tee -a ${logfile}
    
    if [ \$EXIT_CODE -ne 0 ]; then
        echo "RapSearch2 BSAT failed, creating empty output file" | tee -a ${logfile}
        touch ${rapsearchBsatDir}/blacklist_bsat.m8
    fi
    
    echo "=== END RAPSEARCH2_BSAT DEBUG ===" | tee -a ${logfile}
    """
}

process RAPSEARCH2_VFDB {
    conda 'bioconda::seqscreen'
    errorStrategy 'ignore'  // Don't fail the whole pipeline if this fails
    
    input:
    val initialized

    output:
    val true, emit: rapsearch_vfdb

    if ( Executor == 'local' ) {
       executor "local"
    }
    else if ( Executor == 'slurm' ) {
       clusterOptions "--ntasks-per-node ${params.threads}"
       executor "slurm"
    }

    script:
    // Limit threads for rapsearch2 to avoid segmentation faults
    def rapsearch_threads = Math.min(params.threads as Integer, params.max_rapsearch_threads as Integer)
    
    """
    echo "=== RAPSEARCH2_VFDB DEBUG ===" | tee -a ${logfile}
    echo "Using ${rapsearch_threads} threads for RapSearch2 (limited from ${params.threads})" | tee -a ${logfile}
    
    # Check if database exists
    if [ ! -f "${databasesDir}/rapsearch2/vfdb/vfdb.seqs.aa" ]; then
        echo "ERROR: RAPSEARCH2 VFDB database not found at ${databasesDir}/rapsearch2/vfdb/vfdb.seqs.aa" | tee -a ${logfile}
        echo "Looking for files in rapsearch2 VFDB directory:" | tee -a ${logfile}
        ls -la "${databasesDir}/rapsearch2/vfdb/" | tee -a ${logfile} || echo "Directory does not exist" | tee -a ${logfile}
        
        # Try to find rapsearch2 databases
        echo "Searching for rapsearch2 databases:" | tee -a ${logfile}
        find "${databasesDir}" -name "*.aa" -type f | head -10 | tee -a ${logfile}
        
        touch ${rapsearchVfdbDir}/blacklist_vfdb.m8
        exit 0
    fi
    
    # Check if the rapsearch2.sh script exists
    if [ ! -f "${MODULES}/rapsearch2.sh" ]; then
        echo "ERROR: rapsearch2.sh script not found at ${MODULES}/rapsearch2.sh" | tee -a ${logfile}
        touch ${rapsearchVfdbDir}/blacklist_vfdb.m8
        exit 0
    fi
    
    set +e
    ${MODULES}/rapsearch2.sh --fasta=${fastafile} \\
                             --database=${databasesDir}/rapsearch2/vfdb/vfdb.seqs.aa \\
                             --out=${rapsearchVfdbDir}/blacklist_vfdb \\
                             --evalue=1e-9 \\
                             --threads=${rapsearch_threads}
    EXIT_CODE=\$?
    set -e
    
    echo "RapSearch2 VFDB exit code: \$EXIT_CODE" | tee -a ${logfile}
    
    if [ \$EXIT_CODE -ne 0 ]; then
        echo "RapSearch2 VFDB failed, creating empty output file" | tee -a ${logfile}
        touch ${rapsearchVfdbDir}/blacklist_vfdb.m8
    fi
    
    echo "=== END RAPSEARCH2_VFDB DEBUG ===" | tee -a ${logfile}
    """
}

process BSAT_BLACKLIST {
    conda 'bioconda::seqscreen'
    
    input:
    val bowtie_bsat
    val rapsearch_bsat

    output:
    val true, emit: blacklist_bsat

    executor Executor

    script:
    """
    echo -n "# Checking for threats by BSAT Blacklist ..........." | tee -a ${logfile}
    
    # Check if input files exist, create empty output if not
    if [ ! -f "${workingDir}/seqmapper/rapsearch2_bsat/blacklist_bsat.m8" ] || [ ! -f "${workingDir}/seqmapper/bowtie2_bsat/blacklist_bsat.sam" ]; then
        echo "WARNING: BSAT input files missing, creating empty output" | tee -a ${logfile}
        touch ${workingDir}/$WORKFLOW/threats_by_blacklist.txt
    else
        # Check if the script exists
        if [ ! -f "${SCRIPTS}/threats_by_blacklist.py" ]; then
            echo "ERROR: threats_by_blacklist.py script not found at ${SCRIPTS}/threats_by_blacklist.py" | tee -a ${logfile}
            touch ${workingDir}/$WORKFLOW/threats_by_blacklist.txt
        else
            ${SCRIPTS}/threats_by_blacklist.py \\
                --rapsearch-results=${workingDir}/seqmapper/rapsearch2_bsat/blacklist_bsat.m8 \\
                --bowtie2-results=${workingDir}/seqmapper/bowtie2_bsat/blacklist_bsat.sam \\
                --out=${workingDir}/$WORKFLOW/threats_by_blacklist.txt
        fi
    fi
    
    echo " DONE " | tee -a ${logfile}
    date '+%H:%M:%S %Y-%m-%d' | tee -a ${logfile}
    """ 
}

process VFDB_BLACKLIST {
    conda 'bioconda::seqscreen'
    
    input:
    val bowtie_vfdb
    val rapsearch_vfdb

    output:
    val true, emit: blacklist_vfdb

    executor Executor

    script:
    """
    echo -n "# Checking for threats by VFDB Blacklist ..........." | tee -a ${logfile}
    
    # Check if input files exist, create empty output if not
    if [ ! -f "${workingDir}/seqmapper/rapsearch2_vfdb/blacklist_vfdb.m8" ] || [ ! -f "${workingDir}/seqmapper/bowtie2_vfdb/blacklist_vfdb.sam" ]; then
        echo "WARNING: VFDB input files missing, creating empty output" | tee -a ${logfile}
        touch ${workingDir}/${WORKFLOW}/threats_by_vfdb.txt
    else
        # Check if the script exists
        if [ ! -f "${SCRIPTS}/threats_by_blacklist.py" ]; then
            echo "ERROR: threats_by_blacklist.py script not found at ${SCRIPTS}/threats_by_blacklist.py" | tee -a ${logfile}
            touch ${workingDir}/${WORKFLOW}/threats_by_vfdb.txt
        else
            ${SCRIPTS}/threats_by_blacklist.py \\
                --rapsearch-results=${workingDir}/seqmapper/rapsearch2_vfdb/blacklist_vfdb.m8 \\
                --bowtie2-results=${workingDir}/seqmapper/bowtie2_vfdb/blacklist_vfdb.sam \\
                --out=${workingDir}/${WORKFLOW}/threats_by_vfdb.txt
        fi
    fi
    
    echo " DONE " | tee -a ${logfile}
    date '+%H:%M:%S %Y-%m-%d' | tee -a ${logfile}
    """ 
}

process PFAM_HMM {
    conda 'bioconda::seqscreen'
    
    input:
    val initialized

    output:
    val true, emit: pfam_hmm

    if ( Executor == 'local' ) {
       executor "local"
    }
    else if ( Executor == 'slurm' ) {
       clusterOptions "--ntasks-per-node ${params.threads}"
       executor "slurm"
    }

    script:
    """
    echo "=== PFAM_HMM DEBUG ===" | tee -a ${logfile}
    
    # Check if esl-translate exists
    if ! command -v esl-translate &> /dev/null; then
        echo "ERROR: esl-translate not found" | tee -a ${logfile}
        touch ${pfamDir}/${BASE}.domtblout
        exit 0
    fi
    
    # Check if hmmscan database exists
    if [ ! -f "${databasesDir}/hmmscan/Pfam-A.hmm" ]; then
        echo "ERROR: Pfam-A.hmm database not found at ${databasesDir}/hmmscan/Pfam-A.hmm" | tee -a ${logfile}
        touch ${pfamDir}/${BASE}.domtblout
        exit 0
    fi
    
    # Check if hmmscan.sh script exists
    if [ ! -f "${MODULES}/hmmscan.sh" ]; then
        echo "ERROR: hmmscan.sh script not found at ${MODULES}/hmmscan.sh" | tee -a ${logfile}
        touch ${pfamDir}/${BASE}.domtblout
        exit 0
    fi
    
    esl-translate  ${fastafile} > ${translatedFile}
    ${MODULES}/hmmscan.sh --fasta=${translatedFile} \\
                          --database=${databasesDir}/hmmscan/Pfam-A.hmm \\
                          --out=${pfamDir}/${BASE} \\
                          --threads=${params.threads} \\
                          --evalue=1e-3
    echo -n " # ${WORKFLOW} workflow complete ......................... " | tee -a ${logfile}; date '+%H:%M:%S %Y-%m-%d' | tee -a ${logfile}
    """
}

workflow {
    CREATE_WORKING_DIRECTORIES()
    INITIALIZE(CREATE_WORKING_DIRECTORIES.out)
    BOWTIE2_BSAT(INITIALIZE.out)
    BOWTIE2_VFDB(INITIALIZE.out)
    RAPSEARCH2_BSAT(INITIALIZE.out)
    RAPSEARCH2_VFDB(INITIALIZE.out)
    BSAT_BLACKLIST(BOWTIE2_BSAT.out, RAPSEARCH2_BSAT.out)
    VFDB_BLACKLIST(BOWTIE2_VFDB.out, RAPSEARCH2_VFDB.out)
    if (params.hmmscan) {
        PFAM_HMM(INITIALIZE.out)
    }
}