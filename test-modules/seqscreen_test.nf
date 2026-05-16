#!/usr/bin/env nextflow

include { convert_to_nfcore_tuple } from '../subworkflows/local/utils/nf-core-compatibility.nf'
include { SEQSCREEN } from "../modules/local/seqscreen/main.nf"
// note: path of module is relative to the directory containing this file! (./testing/)

// -------------------------
// Parameters
// -------------------------
// note: paths are relative to the directory containing this file
params.fasta = "${projectDir}/../assets/data/other_tools_files/more_files/seqscreen/example.fasta"
params.seqscreen_db = "/home/dbs/SeqScreenDB_23.4/"
params.seqscreen_mode = "fast"

// -------------------------
// Workflow
// -------------------------
workflow {

    fasta = convert_to_nfcore_tuple(params.fasta)
    db = channel.fromPath(params.seqscreen_db)

    SEQSCREEN(fasta, db)
}
    