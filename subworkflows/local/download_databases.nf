#!/usr/bin/env nextflow

include { HOSTILE_FETCH } from '../../modules/nf-core/hostile/fetch/main.nf'
include { CHECKM2_DATABASEDOWNLOAD } from '../../modules/nf-core/checkm2/databasedownload/main'   
include { BAKTA_BAKTADBDOWNLOAD } from '../../modules/nf-core/bakta/baktadbdownload/main' 
include { SINGLEM_DOWNLOAD_DB } from '../../modules/local/singlem/download_db/main.nf'
include { EMU_DOWNLOAD_DB ; EMU_STAGE_DB } from "../../modules/local/emu/download_db/main.nf"
include { LEMUR_DATABASEDOWNLOAD ; LEMUR_STAGE_DB } from "../../modules/local/lemur/download_db/main.nf"
include { SYLPH_DOWNLOAD_DB } from "../../modules/local/sylph/download_db/main.nf"

workflow DOWNLOAD_DBS {

    take:
    analysis_type // string: type of analysis (e.g. 'assembly', 'taxonomic-profiling', 'genome-dynamics')

    hostile_index // string: name of the hostile database file *with extension* (e.g. 'human-t2t-hla-argos985-mycob140.mmi')
    lemur_db_zenodo_id // string: Zenodo ID of the lemur database
    checkm2_db_zenodo_id // string: Zenodo ID of the checkm2 database
    target_taxonomy // string: target taxonomy for the database
    sylph_db_picker // string: path to the CSV file containing database versions


    main:
    // Initialize empty channels for each database type
    ch_hostile_db = channel.empty()
    ch_emu_db     = channel.empty()
    ch_lemur_db   = channel.empty()
    ch_sylph_db   = channel.empty()
    ch_checkm2_db = channel.empty()
    ch_bakta_db   = channel.empty()
    ch_singlem_db = channel.empty()


    // log message: downloading databases for which analysis type
    log.info "Downloading databases for analysis type: ${analysis_type}"


    // ------------------------------------------------
    // pre-processing databases 
    // ------------------------------------------------
    // TODO: need to add conditional for running hostile based on params.sample_environment
    db_name_without_extension = hostile_index.replaceAll('\\.mmi$', '')
    log.info "Fetching hostile index/database: ${db_name_without_extension} for minimap2. Will take > 5 min"
    HOSTILE_FETCH(db_name_without_extension)
    ch_hostile_db = HOSTILE_FETCH.out.reference

    // ------------------------------------------------
    // taxonomic profiling databases 
    // ------------------------------------------------
    if (analysis_type == "taxonomic-profiling") {
        
        if (params.data_type == "16S") {
            // download emu db
            log.info "Downloading EMU database for 16S profiling"
            EMU_DOWNLOAD_DB()
            EMU_STAGE_DB(EMU_DOWNLOAD_DB.out.db_files)
            ch_emu_db = EMU_STAGE_DB.out.emu_db

        } else if (params.taxonomic_profiler == "sylph") {
            log.info "Downloading Sylph database for taxonomic profiling"
            SYLPH_DOWNLOAD_DB(target_taxonomy, sylph_db_picker)
            ch_sylph_db = SYLPH_DOWNLOAD_DB.out.sylph_db
            
        } else {
            // download lemur db
            log.info "Downloading LEMUR database for taxonomic profiling"
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
        log.info "Downloading CheckM2 database for assembly"
        CHECKM2_DATABASEDOWNLOAD(checkm2_db_zenodo_id)
        ch_checkm2_db = CHECKM2_DATABASEDOWNLOAD.out.database
    
        // download bakta db
        log.info "Downloading Bakta database for assembly"
        BAKTA_BAKTADBDOWNLOAD()
        log.warn "If downloading Bakta database which is ~55GB size: it takes ~50 minutes"
        ch_bakta_db = BAKTA_BAKTADBDOWNLOAD.out.db

        // download singlem db
        log.info "Downloading Singlem database for assembly"
        SINGLEM_DOWNLOAD_DB()
        ch_singlem_db = SINGLEM_DOWNLOAD_DB.out.singlem_db
    }

    emit: // emit empty channels if not downloaded
    ch_hostile_db = ch_hostile_db

    // taxonomic profiling databases
    ch_emu_db = ch_emu_db // not used currently ; 
    ch_lemur_db = ch_lemur_db
    ch_sylph_db = ch_sylph_db
    
    // assembly databases
    ch_checkm2_db = ch_checkm2_db
    ch_bakta_db = ch_bakta_db
    ch_singlem_db = ch_singlem_db
}
