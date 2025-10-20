#!/usr/bin/env nextflow

// enable dsl2 syntax
nextflow.enable.dsl = 2

include { convert_to_nfcore_tuple } from '../subworkflows/local/utils/nf-core-compatibility.nf'
include { LEMUR_DATABASEDOWNLOAD ; LEMUR_STAGE_DB } from "../modules/local/lemur/databasedownload/main.nf"
include { LEMUR } from "../modules/local/lemur/main.nf"

// -------------------------
// Parameters
// -------------------------
// Note: Paths are relative to the base directory of the workflow (where nextflow is run from)

params.reads = "${projectDir}/../assets/examples/other_tool_files/lemur/example-data/example.fastq"

// -------------------------
// Workflow
// -------------------------
workflow {
    
    reads_ch = convert_to_nfcore_tuple(params.reads)
    
    LEMUR_DATABASEDOWNLOAD(params.lemur_db_zenodo_id)
    LEMUR_STAGE_DB(LEMUR_DATABASEDOWNLOAD.out.db_files, LEMUR_DATABASEDOWNLOAD.out.refseq_version_bacteria)
    
    LEMUR_STAGE_DB.out.lemur_db.view { x -> "lemur_db: $x"}
    
    // LEMUR(reads_ch)
}