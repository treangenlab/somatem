#!/usr/bin/env nextflow

include { HOSTILE_FETCH } from "../../modules/nf-core/hostile/fetch/main.nf"
include { HOSTILE_CLEAN } from "../../modules/nf-core/hostile/clean/main.nf"
include { convert_to_nfcore_tuple } from "./utils/nf-core-compatibility.nf"
// note: path of module is relative to the directory containing this file! (./testing/)

// -------------------------
// Required Parameters
// -------------------------
// params.hostile_database_dir = "${launchDir}/databases/hostile/reference"

// -------------------------
// Workflow
// Look for/download .mmi (minimap2) database in the database directory and run Hostile
// -------------------------
workflow runHostile {

    take:
    reads // tuple: [ meta, reads] as channels
    database_name // string: name of the database file *with extension* (e.g. 'human-t2t-hla-argos985-mycob140.mmi')

    main:
    // Check if database exists in the default directory
    db_path = "${params.hostile_db}/${database_name}"
    db_exists = Channel.value(file(db_path).exists())

    
    if ( db_exists ) {
        database_input = Channel.value([database_name, params.hostile_db])
    } else { 
        // If database is not in the default directory, fetch it (warning: Takes a while)
        db_name_without_extension = database_name.replaceAll('\\.mmi$', '')
        log.info "Fetching hostile index/database: ${db_name_without_extension} for minimap2. Will take > 5 min"
        HOSTILE_FETCH(db_name_without_extension)
        
        MOVE_DB(HOSTILE_FETCH.out.reference) // move the fetched database to the default directory (shared) for future use
        database_input = MOVE_DB.out.db_tuple

        database_input.map { db_name, db_dir -> 
            log.info "Database: ${db_name} moved to the global database directory: ${db_dir}"
        }
    }
    
    // run Hostile clean
    HOSTILE_CLEAN(reads, database_input)
    
    ch_versions = HOSTILE_CLEAN.out.versions.first() // collect version of hostile
    
    emit:
    dehosted_reads = HOSTILE_CLEAN.out.fastq // de-hosted fastq file
    versions = ch_versions
}
    

// Move downloaded dbs from the work/../reference into db_base_dir/hostile/ dir
process MOVE_DB {
    input:
    tuple val(database_name), path(reference)
    output:
    tuple val(database_name_mmi), path(params.hostile_database_dir), emit: db_tuple
    
    script:
    database_name_mmi = "${database_name}.mmi"

    """    
    cp "${reference}/${database_name_mmi}" "${params.hostile_database_dir}/${database_name_mmi}"
    """
}