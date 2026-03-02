#!/usr/bin/env nextflow

include { HOSTILE_FETCH } from '../../modules/nf-core/hostile/fetch/main.nf'
include { CHECKM2_DATABASEDOWNLOAD } from '../../modules/nf-core/checkm2/databasedownload/main'   
include { BAKTA_BAKTADBDOWNLOAD } from '../../modules/nf-core/bakta/baktadbdownload/main' 
include { SINGLEM_DOWNLOAD_DB } from '../../modules/local/singlem/download_db/main.nf'
include { EMU_DOWNLOAD_DB ; EMU_STAGE_DB } from "../../modules/local/emu/downloaddb/main.nf"
include { LEMUR_DATABASEDOWNLOAD ; LEMUR_STAGE_DB } from "../../modules/local/lemur/databasedownload/main.nf"

workflow DOWNLOAD_DBS {

    take:
    analysis_type // string: type of analysis (e.g. 'assembly', 'taxonomic-profiling', 'genome-dynamics')

    hostile_index // string: name of the hostile database file *with extension* (e.g. 'human-t2t-hla-argos985-mycob140.mmi')
    lemur_db_zenodo_id // string: Zenodo ID of the lemur database
    checkm2_db_zenodo_id // string: Zenodo ID of the checkm2 database


    main:
    // Initialize empty channels for each database type
    ch_hostile_db = Channel.empty()
    ch_emu_db = Channel.empty()
    ch_lemur_db = Channel.empty()
    ch_checkm2_db = Channel.empty()
    ch_bakta_db = Channel.empty()
    ch_singlem_db = Channel.empty()

    // log message: downloading databases for which analysis type
    log.info "Downloading databases for analysis type: ${analysis_type}"


    // ------------------------------------------------
    // pre-processing databases 
    // ------------------------------------------------

    if (params.run_hostile) {
        db_name_without_extension = hostile_index.replaceAll('\\.mmi$', '')

        // TODO: move this log message to the process itself
        log.info "Fetching hostile index/database: ${db_name_without_extension} for minimap2. Will take > 5 min"
        HOSTILE_FETCH(db_name_without_extension)
        ch_hostile_db = HOSTILE_FETCH.out.reference
    } else {
        // skip hostile database download ; ch_hostile_db remains empty
    }


    

    // ------------------------------------------------
    // taxonomic profiling databases 
    // ------------------------------------------------
    if (analysis_type == "taxonomic-profiling") {
        
        if (params.data_type == "16S") {
            // download emu db
            EMU_DOWNLOAD_DB()
            EMU_STAGE_DB(EMU_DOWNLOAD_DB.out.db_files)
            ch_emu_db = EMU_STAGE_DB.out.emu_db

        } else {
            // download lemur db
            LEMUR_DATABASEDOWNLOAD(lemur_db_zenodo_id)
            LEMUR_STAGE_DB(LEMUR_DATABASEDOWNLOAD.out.db_files, LEMUR_DATABASEDOWNLOAD.out.refseq_version_bacteria)
            ch_lemur_db = LEMUR_STAGE_DB.out.lemur_db
        }
    }

    // ------------------------------------------------
    // assembly databases 
    // ------------------------------------------------
    if (analysis_type == "assembly") {
        // download checkm2 database 
        CHECKM2_DATABASEDOWNLOAD(checkm2_db_zenodo_id)
        ch_checkm2_db = CHECKM2_DATABASEDOWNLOAD.out.database
    
        // download bakta db
        BAKTA_BAKTADBDOWNLOAD()
        log.warn "If downloading Bakta database which is ~55GB size: it takes ~50 minutes"
        ch_bakta_db = BAKTA_BAKTADBDOWNLOAD.out.db

        // download singlem db
        SINGLEM_DOWNLOAD_DB()
        ch_singlem_db = SINGLEM_DOWNLOAD_DB.out.singlem_db
    }

    emit: // emit empty channels if not downloaded
    ch_hostile_db = ch_hostile_db

    // taxonomic profiling databases
    ch_emu_db = ch_emu_db // not used currently ; 
    ch_lemur_db = ch_lemur_db
    
    // assembly databases
    ch_checkm2_db = ch_checkm2_db
    ch_bakta_db = ch_bakta_db
    ch_singlem_db = ch_singlem_db
}
