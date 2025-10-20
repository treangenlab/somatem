#!/usr/bin/env nextflow

include { HOSTILE_FETCH } from "../modules/nf-core/hostile/fetch/main.nf"
include { HOSTILE_CLEAN } from "../modules/nf-core/hostile/clean/main.nf"
include { convert_to_nfcore_tuple } from "../subworkflows/local/utils/nf-core-compatibility.nf"
// note: path of module is relative to the directory containing this file! (./testing/)

// -------------------------
// Parameters
// -------------------------
// note: paths are relative to the workflow directory (from where nextflow is run)
params.reads = "${projectDir}/../assets/examples/data/mock9_sub10k.fastq.gz"
params.database_name = 'human-t2t-hla-argos985-mycob140.mmi'
params.database_dir = "${params.db_base_dir}/hostile"

// -------------------------
// Workflow
// -------------------------
workflow {

    reads = convert_to_nfcore_tuple(params.reads)
    // database_name = Channel.of(params.database_name)

    // testing db fetch with storeDir 
    db_name_without_extension = params.database_name.replaceAll('\\.mmi$', '')
    HOSTILE_FETCH(db_name_without_extension)
    
    HOSTILE_CLEAN(reads, HOSTILE_FETCH.out.reference)


    // // if database is in the default directory, use it
    // if ( params.database_name == 'human-t2t-hla-argos985-mycob140.mmi' ) {

    //     database_dir = Channel.of(params.database_dir)
    //     database_tuple = database_name.combine(database_dir)

    //     HOSTILE_CLEAN(reads, database_tuple)

    // } else { 
    //     // if database is not in the default directory, fetch it

    //     printf("fetching database: %s", params.database_name)
    //     HOSTILE_FETCH(database_name)
    //     HOSTILE_CLEAN(reads, HOSTILE_FETCH.out.reference)
    // }
}
    

// tasks
// Need to add ext.args to download only mmi (minimap2) database for HOSTILE_FETCH (with --minimap2 flag)