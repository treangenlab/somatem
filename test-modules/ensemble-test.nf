#!/usr/bin/env nextflow

include { convert_to_nfcore_tuple } from '../subworkflows/local/utils/nf-core-compatibility.nf'
include { SYLPH_PROFILE } from "../modules/nf-core/sylph/profile/main.nf"
include { GANON_CLASSIFY } from "../modules/nf-core/ganon/classify/main.nf"
include { KRAKEN2_KRAKEN2 } from "../modules/nf-core/kraken2/kraken2/main.nf"



// -------------------------
// Parameters
// -------------------------
// Note: Paths are relative to the base directory of the workflow (where nextflow is run from)

// params.reads = "${projectDir}/../assets/examples/other_tool_files/lemur/example-data/example.fastq"
params.reads = './assets/examples/other_tools_files/sylph/o157_reads.fastq.gz'

// -------------------------
// Workflow
// -------------------------
workflow {
    
    reads_ch = convert_to_nfcore_tuple(params.reads)
    
    SYLPH_PROFILE(reads_ch, params.sylph_db)
    GANON_CLASSIFY(reads_ch, params.ganon_db)
    KRAKEN2_KRAKEN2(reads_ch, params.kraken2_db, true, true)
    
    // LEMUR(reads_ch)
}