#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.fasta = ''
params.working = ''
params.bin_dir = null
params.module_dir = null
params.databases = ''
params.threads = 1
params.evalue = 10
params.log = '/dev/null'
params.help = false
params.slurm = false
params.blastn = false
Executor = 'local'

MODULES = "${params.module_dir}"
SCRIPTS = "${params.bin_dir}"

def usage() {
    log.info ""
    log.info "Usage: nextflow run taxonomic_identification.nf --fasta /Path/to/infile.fasta --databases /Path/to/databases --working /Path/to/working_directory [--evalue 10] [--threads 4] [--slurm]  [--log /Path/to/log.txt]"
    log.info "  --fasta     Path to the input NT FASTA file"
    log.info "  --databases Path to the database directory"
    log.info "  --working   Path to the output working directory"
    log.info "  --blastn    Run blastn in addition to blastx"
    log.info "  --evalue    E value cut-off (Default=10)"
    log.info "  --threads   Number of threads to use (Default=1)"
    log.info '  --slurm     Submit modules in this workflow to run on a SLURM grid (Default = run locally)'
    log.info "  --log       Where to write log file (Default=no log)"
    log.info "  -h, --help  Print this help message out"
    log.info ""
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

fastafile = file(params.fasta)
workingDir = file(params.working)
databaseDir = file(params.databases)
logfile = file(params.log)

BASE = fastafile.getName()
THREADS = params.threads
EVALUE = params.evalue
WORKFLOW = "taxonomic_identification"

blastnDir = file("${workingDir}/${WORKFLOW}/blastn")
blastxDir = file("${workingDir}/${WORKFLOW}/blastx")
taxDir = file("${workingDir}/${WORKFLOW}/taxonomic_assignment")
outlierDir = file("${workingDir}/${WORKFLOW}/outlier_detection")

translatedFile = "${workingDir}/initialize/six_frame_translation/${BASE}.translated.fasta"

process CREATE_WORKING_DIRECTORIES {
    conda 'bioconda::seqscreen'

    output:
        val true, emit: working_dir

    script:
    """
    mkdir -p ${workingDir}/${WORKFLOW}/
    if [ -d ${blastnDir} ]; then rm -rf ${blastnDir}; fi;
    if [ -d ${outlierDir} ]; then rm -rf ${outlierDir}; fi;
    if [ -d ${blastxDir} ]; then rm -rf ${blastxDir}; fi;
    if [ -d ${taxDir} ]; then rm -rf ${taxDir}; fi;
    
    mkdir ${blastnDir}
    mkdir ${outlierDir}
    mkdir ${blastxDir}
    mkdir ${taxDir}
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
    echo -n " # Launching sensitive ${WORKFLOW} workflow ...... " | tee -a ${logfile}; date '+%H:%M:%S %Y-%m-%d' | tee -a ${logfile}
    """
}

process BLASTN {
    conda 'bioconda::seqscreen'

    input:
        val initialized

    output:
        val true, emit: blastN

    if ( Executor == 'local' ) {
       executor "local"
    }

    else if ( Executor == 'slurm' ) {
       clusterOptions "--ntasks-per-node ${THREADS}"
       executor "slurm"
    }
    
    script:
    if (params.blastn) {
        """
        echo -n " # Running blastN ...... " | tee -a ${logfile}; date '+%H:%M:%S %Y-%m-%d' | tee -a ${logfile}

        ${MODULES}/blastn.sh --fasta=$fastafile \
                     --database=${databaseDir}/blast/nt/nt \
                 --out=${blastnDir}/${BASE}.nt.btab \
                     --threads=${THREADS} \
                 --evalue=${EVALUE}
        """
    } else {
        """
        """
    }
}

process OUTLIER_DETECTION {
    conda 'bioconda::seqscreen'

    input:
        val blastN

    output:
        val true, emit: outlier

    script:
    if (params.blastn) {
        """
        echo -n " # Running outlier detection ...... " | tee -a ${logfile}; date '+%H:%M:%S %Y-%m-%d' | tee -a ${logfile}
        ${MODULES}/outlier_detection.sh --fasta=$fastafile \
                                        --btab=${blastnDir}/${BASE}.nt.btab \
                                        --out=${outlierDir}/outlier_detection.txt
        """
    } else {
        """
        """
    }
}

process BLASTX {
    conda 'bioconda::seqscreen'

    input:
        val blastN

    output:
        val true, emit: blastX

    if ( Executor == 'local' ) {
       executor "local"
    }

    else if ( Executor == 'slurm' ) {
       clusterOptions "--ntasks-per-node ${THREADS}"
       executor "slurm"
    }

    script:
    """
    echo -n " # Running blastX ...... " | tee -a ${logfile}; date '+%H:%M:%S %Y-%m-%d' | tee -a ${logfile}
    ${MODULES}/blastx.sh --fasta=$fastafile \
    			      --database=${databaseDir}/blast/UNIREF100.mini \
    			      --out=${blastxDir}/${BASE}.ur100 \
    			      --threads=${THREADS} \
    			      --evalue=${EVALUE}
    ln -s ${BASE}.ur100.btab ${blastxDir}/functional_link.ur100.btab
    ln -s ${BASE}.ur100.xml ${blastxDir}/functional_link.ur100.xml
    """
}

process TAXONOMIC_ASSIGNMENT {
    conda 'bioconda::seqscreen'
    
    input:
        val blastN
        val outlier
        val blastX

    output:
        val true, emit: taxonomic_assignment

    executor Executor
    
    script:
    """
    echo -n " # Running sensitive taxonomic assignment ...... " | tee -a ${logfile}; date '+%H:%M:%S %Y-%m-%d' | tee -a ${logfile}
    ${SCRIPTS}/btab_2_tax_report_sensitive.pl --blastx=${blastxDir}/${BASE}.ur100.btab \
    				    --blastn=${blastnDir}/${BASE}.nt.btab_outlier_clean.btab \
    				    --out=${taxDir}/taxonomic_results.txt \
    				    --cutoff=1
    """
}

workflow {
    CREATE_WORKING_DIRECTORIES()
    INITIALIZE(CREATE_WORKING_DIRECTORIES.out)
    BLASTN(INITIALIZE.out)
    OUTLIER_DETECTION(BLASTN.out)
    BLASTX(BLASTN.out)
    TAXONOMIC_ASSIGNMENT(BLASTN.out, OUTLIER_DETECTION.out, BLASTX.out)
}