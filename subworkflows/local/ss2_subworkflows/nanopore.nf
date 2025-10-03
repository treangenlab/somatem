#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.fasta = ""
params.working = ""
params.bin_dir = null
params.module_dir = null
params.databases = ""
params.threads = 1
params.evalue = 0.001
params.splitby = 100000
params.taxonomy_confidence_threshold = 0
params.bitscore_cutoff = false
params.log = '/dev/null'
params.help = false
params.slurm = false
params.ancestral = false
params.includecc = false
params.skip_report = false
params.online = false
params.format = 1
params.keep_html_ont = false
params.filter_taxon = "\"\""
params.keep_taxon = "\"\""
ancestral_flag = ""
include_cc_flag = ""
params.prefix = ''
cutoff_flag = "--cutoff=5"
online_flag = ""
mode="ont"

Executor = 'local'

MODULES = "${params.module_dir}"
SCRIPTS = "${params.bin_dir}"


def usage() {
    log.info ""
    log.info "Usage: nextflow run nanopore.nf --fasta /Path/to/infile.fasta --databases /Path/to/databases --working /Path/to/working_directory [--evalue 10] [--threads 4] [--slurm]  [--log /Path/to/log.txt]"
    log.info "  --fasta     Path to the input NT FASTA file"
    log.info "  --databases Path to the database directory"
    log.info "  --working   Path to the output working directory"
    log.info "  --evalue    E value cut-off (Default=10)"
    log.info "  --threads   Number of threads to use (Default=1)"
    log.info "  --splitby   Number of fasta sequences in each chunk"
    log.info "  --slurm     Submit modules in this workflow to run on a SLURM grid (Default = run locally)"
    log.info "  --ancestral Include all ancestral GO terms" 
    log.info "  --includecc     Include cellular component go terms"
    log.info "  --bitscore_cutoff Tiebreak across all uniprots within this % of the top bitscore"
    log.info "  --prefix   Add argument to beginning of the output files"
    log.info "  --format        Format type: [1] Original, [2] Hits only, [3] FunSoC only, [4] Gene-Centric,[5] Gene-Centric FunSoC Only [Default: 1]"
    log.info "  --skip_report   Skip report generation step and only generate intermediate files" 
    log.info "  --online        Pull reference genomes from NCBI for reference_inference [Needs web access]"
    log.info "  --filter_taxon  Filter comma separated list of taxon"
    log.info "  --keep_taxon    Keep comma separated list of taxon"
    log.info "  --taxonomy_confidence_threshold Confidence threshold for multi-tax ids (Average) [Default: 0.0]"
    log.info "  --keep_html_ont Keep html report in ont mode [Takes additional memory]" 
    log.info "  --version   Verison number for the seqscreen pipeline"
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
rflag = 'ont'
if (params.bitscore_cutoff) {
    cutoff_flag = "--cutoff=${params.bitscore_cutoff}"
    if(rflag) {
        rflag = rflag + ",bitscore:" + params.bitscore_cutoff
    }else{
        rflag = "bitscore:" + params.bitscore_cutoff
    }

}

if (params.evalue != "") {
   if(rflag) {
       rflag = rflag + ",evalue:" + params.evalue
   }
   else {
       rflag = "evalue:" + params.evalue
   }
}


if (params.includecc) {
    include_cc_flag = "--include-cc"
    if(rflag){
        rflag = rflag + ",includecc"
    }else{
        rflag = 'includecc'
    }

}

if (params.ancestral) {
    ancestral_flag = "--ancestral"
    if(rflag) {
      
       rflag = rflag + ",ancestral"
    }
   else {
       rflag = "ancestral"
   }
}

if (params.splitby != ""){
    if(rflag){
        rflag = rflag + ",splitby:" + params.splitby
    }else{
        rflag = "splitby:" + params.splitby
    }
}

if(! rflag)
{
   rflag = "None"
}
 
if (params.online) {
   online_flag = "--online"
}

if (params.taxonomy_confidence_threshold) {
   taxonomy_confidence_threshold_flag = "--taxonomy_confidence_threshold ${params.taxonomy_confidence_threshold}"
}

fastafile = file(params.fasta)
workingDir = file(params.working)
databaseDir = file(params.databases)
logfile = file(params.log)

BASE = fastafile.getName()
WORKFLOW = "taxonomic_identification"
FWORKFLOW="functional_annotation"
RWORKFLOW="report_generation"

blastnDir = file("${workingDir}/${WORKFLOW}/blastn")
blastxDir = file("${workingDir}/${WORKFLOW}/blastx")
taxDir = file("${workingDir}/${WORKFLOW}/taxonomic_assignment")
fxnDir = file("${workingDir}/${FWORKFLOW}/functional_assignments")

process CREATE_WORKING_DIRECTORIES {
    conda 'bioconda::seqscreen'

    output:
        val true, emit: working_dir

    script:
    """
    mkdir -p ${workingDir}/${WORKFLOW}/
    mkdir -p ${workingDir}/${FWORKFLOW}
    mkdir -p ${workingDir}/${RWORKFLOW}
    if [ -d ${blastnDir} ]; then rm -rf ${blastnDir}; fi;
    if [ -d ${blastxDir} ]; then rm -rf ${blastxDir}; fi;
    if [ -d ${taxDir} ]; then rm -rf ${taxDir}; fi;    
    if [ -d ${fxnDir} ]; then rm -rf ${fxnDir}; fi; 
    
    mkdir ${blastnDir}
    mkdir ${blastxDir}
    mkdir ${taxDir}
    mkdir ${fxnDir}
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
    echo -n " # Launching ${WORKFLOW} workflow ...... " | tee -a ${logfile}; date '+%H:%M:%S %Y-%m-%d' | tee -a ${logfile}
    """
}

process DIAMOND {
    conda 'bioconda::seqscreen'

    input:
        val initialized
        path fasta

    output:
        path "${BASE}.ur100.btab", emit: diamond_btab
        path "${BASE}.ur100.xml", emit: diamond_xml

    if ( Executor == 'local' ) {
       executor "local"
       maxForks 1
    }

    else if ( Executor == 'slurm' ) {
       clusterOptions "--ntasks-per-node ${params.threads}"
       executor "slurm"
    }

    script:
    """
    echo -n " # Running diamond on ${fasta} ...... " | tee -a ${logfile}; date '+%H:%M:%S %Y-%m-%d' | tee -a ${logfile}

    ## Run the BLASTx and produce an ASN file
    diamond blastx -q ${fasta} \
        -d "${databaseDir}/diamond/uniref.mini.dmnd" \
        -o "${BASE}" \
        --evalue "${params.evalue}" \
        --threads "${params.threads}" \
        --block-size 4 \
        --index-chunks 3 \
        --salltitles \
        --fast \
        --frameshift 15 \
        --range-culling \
        --culling-overlap 50 \
        --range-cover 50 \
        --min-orf 30 \
        -k 1 \
        --masking seg \
        -f 100 \
        --log
    ${SCRIPTS}/blast_formatter_nano.pl -a ${BASE}.daa -o ${BASE}.ur100 -f 5,"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore ppos qframe score salltitles"
    """
}

process SPLIT_READS {
    conda 'bioconda::seqscreen'

    input:
        path diamond_btab

    output:
        val true, emit: split_reads
    
    executor Executor

    script:
    """
    echo -n " # Running split_reads ...... " | tee -a ${logfile}; date '+%H:%M:%S %Y-%m-%d' | tee -a ${logfile}
      python3 ${SCRIPTS}/split_reads.py --fasta ${fastafile} --blastx ${blastxDir}/${BASE}.ur100.btab
    """
}

process DIAMOND_UM {
    conda 'bioconda::seqscreen'

    input: 
        val split_reads
    
    output:
        path "${BASE}.ur100.btab", emit: diamond_btab_um
        path "${BASE}.ur100.xml", emit: diamond_xml_um
    
    if ( Executor == 'local' ) {
       executor "local"
       maxForks 1
    }

    else if ( Executor == 'slurm' ) {
       clusterOptions "--ntasks-per-node ${params.threads}"
       executor "slurm"
    }

    script:
    if (file("${blastxDir}/${BASE}.unmapped.fasta").size() == 0)
    """ 
    touch ${blastxDir}/${BASE}.ur100.unmapped.btab 
    touch ${blastxDir}/${BASE}.ur100.unmapped.xml
    touch ${BASE}.ur100.btab
    touch ${BASE}.ur100.xml
    """
    else
    """
    echo -n " # Running diamond_unmapped ...... " | tee -a ${logfile}; date '+%H:%M:%S %Y-%m-%d' | tee -a ${logfile}
    ## Run the BLASTx and produce an ASN file
    diamond blastx -q "${blastxDir}/${BASE}.unmapped.fasta" \
        -d "${databaseDir}/diamond_1_2/uniref_1_2.dmnd" \
        -o "${BASE}" \
        --evalue "${params.evalue}" \
        --threads "${params.threads}" \
        --block-size 4 \
        --index-chunks 3 \
        --salltitles \
        --frameshift 15 \
        --range-culling \
        --culling-overlap 50 \
        --range-cover 50 \
        --min-orf 30 \
        --masking seg \
        -k 1 \
        -f 100 \
        --log
        ${SCRIPTS}/blast_formatter_nano.pl -a ${BASE}.daa -o ${BASE}.ur100 -f 5,"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore ppos qframe score salltitles"
    """
}

process SPLIT_READS_UM {
    conda 'bioconda::seqscreen'

    input:
        path diamond_btab_um

    output:
        val true, emit: split_reads_um
    
    executor Executor

    script:
    """
    echo -n " # Splitting unmapped reads ...... " | tee -a ${logfile}; date '+%H:%M:%S %Y-%m-%d' | tee -a ${logfile}
      python3 ${SCRIPTS}/split_reads.py --fasta=${blastxDir}/${BASE}.unmapped.fasta --blastx=${blastxDir}/${BASE}.ur100.unmapped.btab --unmapped
    """
}

process POST_PROCESS_DIAMOND {
    conda 'bioconda::seqscreen'

    input:
        val initialized
        val split_reads_um
    
    output: 
        val true, emit: post_process_diamond
    
    if (Executor == 'local') {                                        
        executor "local"                                              
    }                                                                 
                                                                        
    else if (Executor == "slurm") {                                   
        clusterOptions "--ntasks-per-node ${params.threads}"                   
        executor "slurm"                                              
    }     

    script:
    """
    echo -n " # Post processing diamond ...... " | tee -a ${logfile}; date '+%H:%M:%S %Y-%m-%d' | tee -a ${logfile}
        cat ${blastxDir}/${BASE}.ur100.unmapped.btab >> ${blastxDir}/${BASE}.ur100.btab
        cat ${blastxDir}/${BASE}.ur100.unmapped.xml >> ${blastxDir}/${BASE}.ur100.xml
        ln -s ${BASE}.ur100.btab ${blastxDir}/functional_link.ur100.btab
        ln -s ${BASE}.ur100.xml  ${blastxDir}/functional_link.ur100.xml
    """
}

process TAXONOMIC_ASSIGNMENT {
    conda 'bioconda::seqscreen'

    input:
        val post_process_diamond

    output:
        val true, emit: taxonomic_assignment

    executor Executor
    
    script:
    """
    while [ ! -f ${workingDir}/taxonomic_identification/blastx/${BASE}.ur100.split.btab ]
    do
        echo -n " # Waiting for ${workingDir}/taxonomic_identification/blastx/${BASE}.ur100.split.btab ...... " | tee -a ${logfile}; date '+%H:%M:%S %Y-%m-%d' | tee -a ${logfile}
        sleep 5
    done
    
    echo -n " # Launching ${WORKFLOW} workflow ...... " | tee -a ${logfile}; date '+%H:%M:%S %Y-%m-%d' | tee -a ${logfile}
    ${SCRIPTS}/btab_2_tax_report_sensitive.pl --blastx=${blastxDir}/${BASE}.ur100.split.btab \
                                              --blastn=${blastnDir}/${BASE}.nt.btab_outlier_clean.btab \
    				                          --out=${taxDir}/taxonomic_results.txt \
    				                          $cutoff_flag
    python3 ${SCRIPTS}/nano_tax.py --blastx=${blastxDir}/${BASE}.ur100.split.btab \
                                   --database=${databaseDir} \
                                   --output=${taxDir}
    echo -n " # ${WORKFLOW} workflow complete ....... " | tee -a ${logfile}; date '+%H:%M:%S %Y-%m-%d' | tee -a ${logfile}
    """
}

process FUNCTIONAL_ASSIGNMENTS {
    conda 'bioconda::seqscreen'

    input:
        val initialized
        val taxonomic_assignment

    output:
        val true, emit: functional_assignments

    executor Executor

    script:
    """
    echo -n " # Launching ${FWORKFLOW} workflow ...... " | tee -a ${logfile}; date '+%H:%M:%S %Y-%m-%d' | tee -a ${logfile}
    while [ ! -f ${workingDir}/taxonomic_identification/blastx/${BASE}.ur100.split.btab ]
    do
      sleep 5
    done

    python3 ${SCRIPTS}/functional_report_generation.py --fasta ${workingDir}/taxonomic_identification/blastx/${BASE}.split.fasta \
                                               --urbtab ${workingDir}/taxonomic_identification/blastx/${BASE}.ur100.split.btab \
                                               --go ${databaseDir}/go/go_network.txt \
                                               --out ${workingDir}/${FWORKFLOW}/functional_assignments/functional_results.txt \
                                               --annotation ${databaseDir}/annotation_scores.pck \
                                               ${ancestral_flag} \
                                               ${include_cc_flag} \
                                               ${cutoff_flag} > ${workingDir}/${FWORKFLOW}/functional_assignments.log 2>&1
    echo -n " # ${FWORKFLOW} workflow complete ....... " | tee -a ${logfile}; date '+%H:%M:%S %Y-%m-%d' | tee -a ${logfile}
    """
}

process SEQSCREEN_TSV_REPORT {
    conda 'bioconda::seqscreen'

    input:
        val initialized
        val functional_assignments

    output:
        val true, emit: seqscreen_tsv_report

    executor Executor

    executor Executor
    
    script:
    if (params.skip_report)
    """
    """
    else 
    """
    echo -n " # Launching ${RWORKFLOW} workflow ...... " | tee -a ${logfile}; date '+%H:%M:%S %Y-%m-%d' | tee -a ${logfile}
    ${SCRIPTS}/seqscreen_tsv_report.py \
                                    --taxonomy=${workingDir}/taxonomic_identification/taxonomic_assignment/taxonomic_results.txt \
                                    --functional=${workingDir}/functional_annotation/functional_assignments/functional_results.txt \
                                    --taxlookup=${databaseDir}/taxonomy/taxa_lookup.txt \
                                    --funsocs=${databaseDir}/funsocs \
                                    --fasta=${workingDir}/taxonomic_identification/blastx/${BASE}.split.fasta \
                                    --mode=${mode} \
                                    --out=${workingDir}/${RWORKFLOW}/${params.prefix}seqscreen_report.tsv \
                                    --parenttax=${databaseDir}/tax_to_parent.pck \
                                    --mergedtax=${databaseDir}/merged_taxa.pck
    """
}

process SEQSCREEN_HTML_REPORT {
    conda 'bioconda::seqscreen'

    input:
        val seqscreen_tsv_report

    output:
        val true, emit: seqscreen_html_report

    executor Executor

    script:
    if (params.keep_html_ont && !params.skip_report)
    """
    echo -n " # Launching html report workflow ...... " | tee -a ${logfile}; date '+%H:%M:%S %Y-%m-%d' | tee -a ${logfile}
    ${MODULES}/html_report_generation.sh \
        --report=${workingDir}/${RWORKFLOW}/${params.prefix}seqscreen_report.tsv \
        --fasta=${workingDir}/taxonomic_identification/blastx/${BASE}.split.fasta \
        --blastx=${workingDir}/taxonomic_identification/blastx/functional_link.ur100.xml \
        --version=${params.version} \
        --mode=${mode} \
        --rflag=${rflag} \
        --gonetwork=${databaseDir}/go/go_network.txt \
        --out=${workingDir}/${RWORKFLOW}/${params.prefix}seqscreen_html_report/
    echo -n " # ${RWORKFLOW} workflow complete ....... " | tee -a ${logfile}; date '+%H:%M:%S %Y-%m-%d' | tee -a ${logfile}
    """
    else
    """
    """
}

process INFERENCE {
    conda 'bioconda::seqscreen'

    input:
        val seqscreen_html_report

    output:
        val true, emit: inference

    executor Executor

    script:
    """
    echo -n "# Launching inference workflow ...... " | tee -a ${logfile}; date '+%H:%M:%S %Y-%m-%d' | tee -a ${logfile}
    python3 ${SCRIPTS}/reference_inference.py --fasta=${fastafile} \
                                              --output=${workingDir} \
                                              --working=${taxDir}/inference_working \
                                              --databases=${databaseDir} \
                                              --threads=${params.threads} \
                                              ${online_flag}
    echo -n "# inference complete ......." | tee -a ${logfile}; date '+%H:%M:%S %Y-%m-%d' | tee -a ${logfile}
    """
}

process FORMAT {
    conda 'bioconda::seqscreen'

    input:
        val inference

    output:
        val true, emit: format

    executor  Executor

    script:
    if (params.skip_report)
    """
    """
    else
    """
    echo -n " # Launching format workflow ...... " | tee -a ${logfile}; date '+%H:%M:%S %Y-%m-%d' | tee -a ${logfile}
    python3 ${SCRIPTS}/format.py --report $workingDir/report_generation \
                                    --format=${params.format} \
                                    --mode ${mode} \
                                    --databases=${databaseDir} \
                                    --taxonomy_confidence_threshold=${params.taxonomy_confidence_threshold} \
                                    --filter-taxon ${params.filter_taxon} \
                                    --keep-taxon ${params.keep_taxon}
    echo -n " # format complete ....... " | tee -a ${logfile}; date '+%H:%M:%S %Y-%m-%d' | tee -a ${logfile}
    """
}

process CLEANUP {
    conda 'bioconda::seqscreen'
    
    input:
        val format
    
    output:
        val true, emit: cleanup

    
    script:
    """
    echo -n " # Launching cleanup workflow ...... " | tee -a ${logfile}; date '+%H:%M:%S %Y-%m-%d' | tee -a ${logfile}
    rm -r ${workingDir}/taxonomic_identification/blastn
    rm -r "${blastxDir}/${BASE}.unmapped.fasta.unmapped.fasta"
    echo -n " # cleanup complete ....... " | tee -a ${logfile}; date '+%H:%M:%S %Y-%m-%d' | tee -a ${logfile}
    """
}

workflow {
    CREATE_WORKING_DIRECTORIES()
    INITIALIZE(CREATE_WORKING_DIRECTORIES.out)
    fasta_ch = Channel.fromPath(fastafile).splitFasta(by: params.splitby, file: true)
    DIAMOND(INITIALIZE.out, fasta_ch)

    diamond_btab_collect = DIAMOND.out.diamond_btab.collectFile(name: "${blastxDir}/${BASE}.ur100.btab")
    DIAMOND.out.diamond_xml.collectFile(name: "${blastxDir}/${BASE}.ur100.xml")

    SPLIT_READS(DIAMOND.out.diamond_btab)
    DIAMOND_UM(SPLIT_READS.out)

    diamond_btab_um_collect = DIAMOND_UM.out.diamond_btab_um.collectFile(name: "${blastxDir}/${BASE}.ur100.unmapped.btab")
    diamond_xml_um_collect = DIAMOND_UM.out.diamond_xml_um.collectFile(name: "${blastxDir}/${BASE}.ur100.unmapped.xml")

    SPLIT_READS_UM(DIAMOND_UM.out.diamond_btab_um)
    POST_PROCESS_DIAMOND(INITIALIZE.out, SPLIT_READS_UM.out)

    TAXONOMIC_ASSIGNMENT(POST_PROCESS_DIAMOND.out)
    FUNCTIONAL_ASSIGNMENTS(INITIALIZE.out, TAXONOMIC_ASSIGNMENT.out)
    SEQSCREEN_TSV_REPORT(INITIALIZE.out, FUNCTIONAL_ASSIGNMENTS.out)
    SEQSCREEN_HTML_REPORT(SEQSCREEN_TSV_REPORT.out)
    INFERENCE(SEQSCREEN_HTML_REPORT.out)
    FORMAT(INFERENCE.out)
    CLEANUP(FORMAT.out)
}