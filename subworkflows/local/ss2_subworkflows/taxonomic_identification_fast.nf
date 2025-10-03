#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.fasta = ''
params.working = ''
params.bin_dir = null
params.module_dir = null
params.databases = ''
params.threads = 1
params.evalue = 10
params.taxlimit = 25
params.splitby = 100000
params.log = '/dev/null'
params.help = false
params.slurm = false
Executor = 'local'

MODULES = "${params.module_dir}"
SCRIPTS = "${params.bin_dir}"

def usage() {
    log.info ""
    log.info "Usage: nextflow run taxonomic_identification.nf --fasta /Path/to/infile.fasta --databases /Path/to/databases --working /Path/to/working_directory [--evalue 10] [--threads 4] [--slurm]  [--log /Path/to/log.txt]"
    log.info "  --fasta     Path to the input NT FASTA file"
    log.info "  --databases Path to the database directory"
    log.info "  --working   Path to the output working directory"
    log.info "  --evalue    E value cut-off (Default=10)"
    log.info "  --threads   Number of threads to use (Default=1)"
    log.info "  --splitby   Number of fasta sequences in each chunk"
    log.info "  --taxlimit  Maximum number of taxIDs to output for a single input"
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
    if [ -d ${blastxDir} ]; then rm -rf ${blastxDir}; fi;
    if [ -d ${taxDir} ]; then rm -rf ${taxDir}; fi;
    
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
    echo -n " # Launching fast ${WORKFLOW} workflow ...... " | tee -a ${logfile}; date '+%H:%M:%S %Y-%m-%d' | tee -a ${logfile}
    """
}

process DIAMOND {
    conda 'bioconda::seqscreen'
    
    // FIXED: Add memory and resource directives
    memory '240 GB'
    cpus { THREADS }
    time '4h'

    input:
        val initialized
        path fasta

    output:
        path "${BASE}.ur100.btab", emit: diamond_btab
        path "${BASE}.ur100.xml",  emit: diamond_xml

    if ( Executor == 'local' ) {
       executor "local"
       maxForks 1
    }
    else if ( Executor == 'slurm' ) {
       clusterOptions "--ntasks-per-node ${THREADS} --mem=16G --time=04:00:00"
       executor "slurm"
    }

    script:
    // FIXED: Reduce memory-intensive parameters
    def reduced_threads = Math.min(THREADS as Integer, 32) // Cap threads to reduce memory usage
    def block_size = 2.0  // Reduced from 200 to 2.0 GB
    def index_chunks = 4  // Increased from 1 to 4 to reduce memory per chunk
    """
    echo -n " # Running diamond on ${fasta} ...... " | tee -a ${logfile}; date '+%H:%M:%S %Y-%m-%d' | tee -a ${logfile}

    diamond blastx -q ${fasta} \
        -d "${databaseDir}/diamond/uniref.mini.dmnd" \
        -o ${BASE} \
        --evalue "${EVALUE}" \
        --threads "${reduced_threads}" \
        --block-size ${block_size} \
        --index-chunks ${index_chunks} \
        --salltitles \
        --more-sensitive \
        --min-orf 10 \
        --masking 0 \
        --top 5 \
        -f 100 \
    && ${SCRIPTS}/blast_formatter_fast.pl -a ${BASE}.daa -o ${BASE}.ur100 -f 5,"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore ppos qframe score salltitles"
    """
}

process POST_PROCESS_DIAMOND {
    conda 'bioconda::seqscreen'

    input: 
        path diamond_btab
        path diamond_xml

    output: 
        val true, emit: post_process_diamond

    if (Executor == 'local') {                                        
        executor "local"                                              
    }                                                                 
    else if (Executor == "slurm") {                                   
        clusterOptions "--ntasks-per-node ${THREADS}"                   
        executor "slurm"                                              
    }   

    script: 
    """
    echo -n " # Running post process diamond ...... " | tee -a ${logfile}; date '+%H:%M:%S %Y-%m-%d' | tee -a ${logfile}

    # Force, no error if already exists; link to the *actual* staged files
    ln -sfn ${diamond_btab} ${blastxDir}/functional_link.ur100.btab
    ln -sfn ${diamond_xml}  ${blastxDir}/functional_link.ur100.xml
    """
}

process CENTRIFUGE {
    conda 'bioconda::seqscreen'

    input:
        val initialized
    
    output:
        path "${BASE}.centrifuge", emit: centrifuge_tsv   // CHANGED: produce a local file Nextflow can stage

    if (Executor == 'local') {                                        
        executor "local"                                              
    }                                                                 
    else if (Executor == "slurm") {                                   
        clusterOptions "--ntasks-per-node ${THREADS}"                   
        executor "slurm"                                              
    }                                                                 

    script:
    """
    echo -n " # Running centrifuge ...... " | tee -a ${logfile}; date '+%H:%M:%S %Y-%m-%d' | tee -a ${logfile}
    ${MODULES}/centrifuge.sh --fasta=${fastafile} \
                    --database=${databaseDir}/centrifuge/abv \
                    --out=${BASE}.centrifuge \
                    --threads=${params.threads}

    # Keep your original side-effect file in the blastxDir as well
    mkdir -p ${blastxDir}
    cp -f ${BASE}.centrifuge ${blastxDir}/${BASE}.centrifuge
    """    
}

process MERGE_CENTRIFUGE_DIAMOND {
    conda 'bioconda::seqscreen'
    
    input:
        path diamond_btab                          // CHANGED: actual file
        path centrifuge_tsv                        // CHANGED: actual file

    output:
        val true, emit: merged                     // unchanged; file is still written to ${taxDir}

    executor Executor                                                 

    script:                                                        
    """          
    echo -n " # Merging centrifuge and diamond ...... " | tee -a ${logfile}; date '+%H:%M:%S %Y-%m-%d' | tee -a ${logfile}
    mkdir -p ${taxDir}
    python3 ${SCRIPTS}/consolidatedtax_2_taxreport_fast.py \
        --diamond=${diamond_btab} \
        --centrifuge=${centrifuge_tsv} \
        --out=${taxDir}/taxonomic_results.txt \
        --cutoff=1 \
        --taxlimit=${params.taxlimit}
    """     
}

workflow {
    CREATE_WORKING_DIRECTORIES()
    INITIALIZE(CREATE_WORKING_DIRECTORIES.out)

    fasta_ch = Channel.fromPath(fastafile).splitFasta(by: params.splitby, file: true)

    DIAMOND(INITIALIZE.out, fasta_ch)

    // collect the btab and xml into single files in blastxDir
    diamond_collect_ch = DIAMOND.out.diamond_btab.collectFile(name: "${blastxDir}/${BASE}.ur100.btab")
    diamond_xml_ch     = DIAMOND.out.diamond_xml.collectFile(name: "${blastxDir}/${BASE}.ur100.xml")

    // Change to channels
    POST_PROCESS_DIAMOND(diamond_collect_ch, diamond_xml_ch)

    CENTRIFUGE(INITIALIZE.out)

    // CHANGED: feed MERGE the actual files so it starts only when they exist
    MERGE_CENTRIFUGE_DIAMOND(diamond_collect_ch, CENTRIFUGE.out.centrifuge_tsv)
}