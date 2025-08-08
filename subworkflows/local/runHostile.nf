#!/usr/bin/env nextflow

include { HOSTILE_FETCH } from "../../modules/nf-core/hostile/fetch/main.nf"
include { HOSTILE_CLEAN } from "../../modules/nf-core/hostile/clean/main.nf"
include { convert_to_nfcore_tuple } from "./utils/nf-core-compatibility.nf"
// note: path of module is relative to the directory containing this file! (./testing/)

// -------------------------
// Parameters
// -------------------------
params.database_dir = "${projectDir}/../../databases/hostile/reference"

// -------------------------
// Workflow
// -------------------------
workflow runHostile {

    take:
    reads // tuple: [ meta, reads] as channels
    database_name // string: name of the database file *with extension* (e.g. 'human-t2t-hla-argos985-mycob140.mmi')

    main:
    // if database is in the default directory, use it
    if ( database_name == 'human-t2t-hla-argos985-mycob140.mmi' ) {
        database_input = Channel.value([database_name, params.database_dir])

    } else { 
        // If database is not in the default directory, fetch it (warning: Takes a while)
        log.info "Fetching database: ${database_name}. Will take > 5 min"
        database_name_ch = Channel.value(database_name)
        HOSTILE_FETCH(database_name_ch)
        
        // Make sure this is a value channel if it should be reused
        database_input = HOSTILE_FETCH.out.reference.first()
    }
    
    HOSTILE_CLEAN(reads, database_input)
    
    // note: need to copy the fetched database to the default directory
    // and implement dynamic matching of the database name to the fetched database
    

    emit:
    HOSTILE_CLEAN.out.fastq // de-hosted fastq file
    // note: this is first draft, Need to add versions etc for nf-core compatibility
}
    

// tasks
// Need to add ext.args to download only mmi (minimap2) database for HOSTILE_FETCH (with --minimap2 flag)