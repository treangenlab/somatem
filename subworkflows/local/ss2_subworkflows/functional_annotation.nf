#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.fasta = ''
params.working = ''
params.bin_dir = null
params.module_dir = null
params.databases = ''
params.threads = 1
params.log = '/dev/null'
params.help = false
params.slurm = false
params.ancestral = false
params.sensitive = false
params.bitscore_cutoff = 5
params.evalue = ""
ancestral_flag = ""
include_cc_flag = ""
params.includecc = false

Executor = 'local'

def usage() {
    log.info ''
    log.info 'Usage: nextflow run functional_annotation.nf --fasta /Path/to/infile.fasta --working /Path/to/working_directory --databases /Path/to/databases [--evalue 10] [--threads 4] [--slurm] [--log=/Path/to/run.log]'
    log.info '  --fasta     Path to the input NT FASTA file'
    log.info '  --working   Path to the output working directory'
    log.info "  --databases Path to the SEQSCREEN database directory"
    log.info "  --evalue    E value cut-off (Default=30)"
    log.info "  --threads   Number of threads to use (Default=1)"
    log.info '  --slurm     Submit modules in this workflow to run on a SLURM grid (Default = run locally)'
    log.info "  --ancestral Include all ancestral GO terms" 
    log.info "  --includecc     Include cellular component go terms"
    log.info "  --bitscore_cutoff Tiebreak across all uniprots within this % of the top bitscore"
    log.info "  --sensitive Run in sensitive mode"
    log.info "  --log       Where to write the log file to (Default=No log)"
    log.info '  --help      Print this help message out'
    log.info ''
    exit 1
}

if (params.help) {
    usage()
}

if (params.evalue) {
    EVALUE = params.evalue
} else {
    if (params.sensitive) {
        EVALUE = 10
    } else {
        EVALUE = 30
    }
}

if (params.slurm) {
    Executor = 'slurm'
}

if (!params.fasta) {
    log.error "Missing argument for --fasta option"
    usage()
}

if (params.includecc) {
    include_cc_flag = "--include-cc"
}

if (!params.working) {
    log.error "Missing argument for --working option"
    usage()
}
if (params.ancestral) {
    ancestral_flag = "--ancestral"
}

if (!params.databases) {
    log.error "Missing argument for --databases option"
    usage()
}

fastafile = file(params.fasta)
workingDir = file(params.working)
databasesDir = file(params.databases)
logfile = file(params.log)

BASE=fastafile.getName()
MODULES = "${params.module_dir}"
SCRIPTS = "${params.bin_dir}"
THREADS=params.threads
WORKFLOW="functional_annotation"

mummerDir = file("${workingDir}/${WORKFLOW}/mummer")
megaresDir = file("${workingDir}/${WORKFLOW}/megares")
fxnDir = file("${workingDir}/${WORKFLOW}/functional_assignments")

translatedFile = "${workingDir}/initialize/six_frame_translation/${BASE}.translated.fasta"

process CREATE_WORKING_DIRECTORIES {
    conda 'bioconda::seqscreen'

    output:
        val true, emit: working_dir

    script:
    if (params.sensitive) {
        """
        mkdir -p ${workingDir}/${WORKFLOW}
        if [ -d ${mummerDir} ]; then rm -rf ${mummerDir}; fi;
        if [ -d ${megaresDir} ]; then rm -rf ${megaresDir}; fi;
        if [ -d ${fxnDir} ]; then rm -rf ${fxnDir}; fi;
        mkdir -p ${fxnDir}
        mkdir -p ${mummerDir}
        mkdir -p ${megaresDir}
        """
    } else {
        """
        mkdir -p ${workingDir}/${WORKFLOW}
        if [ -d ${fxnDir} ]; then rm -rf ${fxnDir}; fi;
        mkdir -p ${fxnDir}
        """
    }
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
    echo -n " # Launching ${WORKFLOW} workflow ...... " | tee -a ${logfile}; date '+%H:%M:%S %Y-%m-%d' | tee -a ${logfile}
    """
}

process MUMMER {
    conda 'bioconda::seqscreen'

    input:
        val initialized

    output:
        val true, emit: mummer

    if ( Executor == 'local' ) {
       executor "local"
    }

    else if ( Executor == 'slurm' ) {
       clusterOptions "--ntasks-per-node ${THREADS}"
       executor "slurm"
    }
    
    script:
    if (params.sensitive) {
        """
        echo -n " # Running mummer ...... " | tee -a ${logfile}; date '+%H:%M:%S %Y-%m-%d' | tee -a ${logfile}
        ${MODULES}/mummer.sh --fasta=${fastafile} \
                             --database=${databasesDir}/rebase/rebase.fna \
                             --out=${mummerDir}/${BASE}.mummer_re.txt
        """
    } else {
        """
        echo -n " # Skipping mummer (not in sensitive mode) ...... " | tee -a ${logfile}; date '+%H:%M:%S %Y-%m-%d' | tee -a ${logfile}
        """
    }
}

process MEGARES {
    conda 'bioconda::seqscreen'

    input:
        val initialized

    output:
        val true, emit: megares

    if ( Executor == 'local' ) {
       executor "local"
    }

    else if ( Executor == 'slurm' ) {
       clusterOptions "--ntasks-per-node ${THREADS}"
       executor "slurm"
    }

    script:
    if (params.sensitive) {
        """
        echo -n " # Running megares ...... " | tee -a ${logfile}; date '+%H:%M:%S %Y-%m-%d' | tee -a ${logfile}
        ${MODULES}/blastn.sh --fasta=${fastafile} \
                             --database=${databasesDir}/megares/megares_full \
                             --out=${megaresDir}/${BASE}.megares \
                             --threads=${THREADS} \
                             --evalue=${EVALUE}
        """
    } else { 
        """
        echo -n " # Skipping megares (not in sensitive mode) ...... " | tee -a ${logfile}; date '+%H:%M:%S %Y-%m-%d' | tee -a ${logfile}
        """ 
    }
}

process FUNCTIONAL_ASSIGNMENTS {
    conda 'bioconda::seqscreen'
    
    input:
        val initialized
        val mummer
        val megares

    output:
        val true, emit: functional_assignments

    executor Executor
    
    script:
    """
    echo -n " # Running functional assignments ....... " | tee -a ${logfile}; date '+%H:%M:%S %Y-%m-%d' | tee -a ${logfile}
    
    # Check if required input files exist before running
    URBTAB_FILE="${workingDir}/taxonomic_identification/blastx/functional_link.ur100.btab"
    GO_FILE="${databasesDir}/go/go_network.txt"
    ANNOTATION_FILE="${databasesDir}/annotation_scores.pck"
    
    echo "Checking input files..." | tee -a ${logfile}
    
    if [ ! -f "\$URBTAB_FILE" ]; then
        echo "ERROR: URBTAB file not found: \$URBTAB_FILE" | tee -a ${logfile}
        echo "Available files in taxonomic_identification/blastx/:" | tee -a ${logfile}
        ls -la "${workingDir}/taxonomic_identification/blastx/" 2>/dev/null | tee -a ${logfile} || echo "Directory does not exist" | tee -a ${logfile}
        exit 1
    fi
    
    if [ ! -f "\$GO_FILE" ]; then
        echo "ERROR: GO network file not found: \$GO_FILE" | tee -a ${logfile}
        echo "Available files in go directory:" | tee -a ${logfile}
        ls -la "${databasesDir}/go/" 2>/dev/null | tee -a ${logfile} || echo "Directory does not exist" | tee -a ${logfile}
        exit 1
    fi
    
    if [ ! -f "\$ANNOTATION_FILE" ]; then
        echo "ERROR: Annotation file not found: \$ANNOTATION_FILE" | tee -a ${logfile}
        echo "Available files in database root:" | tee -a ${logfile}
        ls -la "${databasesDir}/" | grep -E "annotation|pck" | tee -a ${logfile}
        exit 1
    fi
    
    echo "All input files found, proceeding with functional assignments..." | tee -a ${logfile}
    
    # Ensure output directory exists
    mkdir -p ${workingDir}/${WORKFLOW}/functional_assignments
    
    # Run the functional report generation
    ${SCRIPTS}/functional_report_generation.py --fasta=${fastafile} \
                                               --urbtab=\$URBTAB_FILE \
                                               --go=\$GO_FILE \
                                               --out=${workingDir}/${WORKFLOW}/functional_assignments/functional_results.txt \
                                               --annotation=\$ANNOTATION_FILE \
                                               ${ancestral_flag} \
                                               ${include_cc_flag} \
                                               --cutoff=${params.bitscore_cutoff} > ${workingDir}/${WORKFLOW}/functional_assignments.log 2>&1
    
    # Check if the command succeeded
    if [ \$? -eq 0 ]; then
        echo -n " # ${WORKFLOW} workflow complete ....... " | tee -a ${logfile}; date '+%H:%M:%S %Y-%m-%d' | tee -a ${logfile}
    else
        echo "ERROR: functional_report_generation.py failed. Check log file: ${workingDir}/${WORKFLOW}/functional_assignments.log" | tee -a ${logfile}
        echo "Last 20 lines of error log:" | tee -a ${logfile}
        tail -20 ${workingDir}/${WORKFLOW}/functional_assignments.log | tee -a ${logfile}
        exit 1
    fi
    """
}

workflow {
    CREATE_WORKING_DIRECTORIES()
    INITIALIZE(CREATE_WORKING_DIRECTORIES.out)
    MUMMER(INITIALIZE.out)
    MEGARES(INITIALIZE.out)
    FUNCTIONAL_ASSIGNMENTS(INITIALIZE.out, MUMMER.out, MEGARES.out)
}