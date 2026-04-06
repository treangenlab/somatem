#!/usr/bin/env nextflow
// this script works when moved into the somatem/ directory and executed with `nextflow run sylph-test.nf`

include { SYLPH_PROFILE } from "./modules/nf-core/sylph/profile/main.nf"
include { convert_to_nfcore_tuple } from './subworkflows/local/utils/nf-core-compatibility.nf'


// -------------------------
// Parameters
// -------------------------
params.reads = './assets/examples/other_tools_files/sylph/o157_reads.fastq.gz'

// params.reference = './databases/v0.3-c1000-gtdb-r214.syldb' // small database
params.threads = 4

// -------------------------
// Workflow
// -------------------------
workflow {

    reads_ch = convert_to_nfcore_tuple(params.reads)

    SYLPH_PROFILE(reads_ch, params.sylph_db)
    // reads.view {read -> "read: ${read.baseName}"}
    // reference.view {ref -> "reference: ${ref.baseName}"}
}
    