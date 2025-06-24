#!/usr/bin/env nextflow

// enable dsl2 syntax
nextflow.enable.dsl = 2

include { lemur } from "./modules/local/lemur/main.nf"

// -------------------------
// Parameters
// -------------------------
params.reads = "./examples/lemur/example-data/example.fastq"

params.database_dir = "./examples/lemur/example-db"
params.taxonomy = "./examples/lemur/example-db/taxonomy.tsv"
params.rank = "species"

// params.output_dir = "./examples/lemur/example-output"


// -------------------------
// Workflow
// -------------------------
workflow {
    
    reads = Channel.fromPath(params.reads)
    database_dir = Channel.fromPath(params.database_dir)
    taxonomy = Channel.fromPath(params.taxonomy)
    rank = Channel.of(params.rank)
    // output_dir = Channel.fromPath(params.output_dir)

    lemur(reads, database_dir, taxonomy, rank)
}