#!/usr/bin/env nextflow

include { DEACON_FETCH } from '../../modules/local/deacon/fetch/main.nf'
include { CHECKM2_DATABASEDOWNLOAD } from '../../modules/nf-core/checkm2/databasedownload/main'   
include { BAKTA_BAKTADBDOWNLOAD } from '../../modules/nf-core/bakta/baktadbdownload/main' 
include { KRAKEN2_STANDARD8_DOWNLOAD_DB; KRAKEN2_STANDARD8_DOWNLOAD_DB as KRAKEN2_TAXONOMY_DOWNLOAD_DB } from '../../modules/local/kraken2/download_db/main.nf'
include { SINGLEM_DOWNLOAD_DB; SINGLEM_DOWNLOAD_DB as SINGLEM_TAXONOMY_DOWNLOAD_DB } from '../../modules/local/singlem/download_db/main.nf'
include { EMU_DOWNLOAD_DB ; EMU_STAGE_DB } from "../../modules/local/emu/download_db/main.nf"
include { LEMUR_DATABASEDOWNLOAD ; LEMUR_STAGE_DB } from "../../modules/local/lemur/download_db/main.nf"
include { SYLPH_DOWNLOAD_DB } from '../../modules/local/sylph/download_db/main.nf'

workflow DOWNLOAD_DBS {

    take:
    analysis_type // string: type of analysis (e.g. 'assembly', 'taxonomic-profiling', 'genome-dynamics')

    deacon_index // string: prebuilt Deacon index name (panhuman-1 or panmouse-1)
    lemur_db_zenodo_id // string: Zenodo ID of the lemur database
    checkm2_db_zenodo_id // string: Zenodo ID of the checkm2 database


    main:
    // Initialize empty channels for each database type
    ch_deacon_index = channel.empty()
    ch_emu_db = channel.empty()
    ch_lemur_db = channel.empty()
    ch_sylph_db = channel.empty()
    ch_checkm2_db = channel.empty()
    ch_bakta_db = channel.empty()
    ch_kraken2_db = channel.empty()
    ch_singlem_db = channel.empty()
    ch_taxonomy_db = channel.empty()

    // log message: downloading databases for which analysis type
    log.info "Downloading databases for analysis type: ${analysis_type}"


    // ------------------------------------------------
    // pre-processing databases 
    // ------------------------------------------------

    if (params.run_deacon && analysis_type != "isolate-analysis" && analysis_type != "seqscreen") {
        if (file(deacon_index).exists()) {
            ch_deacon_index = channel.value(file(deacon_index))
        } else {
            log.info "Fetching Deacon host-depletion index: ${deacon_index}"
            DEACON_FETCH(deacon_index)
            ch_deacon_index = DEACON_FETCH.out.index
        }
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

        } else if (params.taxonomic_profiler == "sylph") {
            if (file(params.sylph_profile_db).exists()) {
                ch_sylph_db = channel.value(file(params.sylph_profile_db))
            } else if (params.sylph_db || params.sylph_taxonomic_group == 'all') {
                error("Sylph database not found: ${params.sylph_profile_db}. Supply an existing database with --sylph_db.")
            } else {
                def sylph_db_name = file(params.sylph_profile_db).name
                def sylph_db_url = "${params.sylph_db_url_base}/${sylph_db_name}"
                SYLPH_DOWNLOAD_DB(sylph_db_url, sylph_db_name)
                ch_sylph_db = SYLPH_DOWNLOAD_DB.out.db
            }
            ch_taxonomy_db = ch_sylph_db
        } else if (params.taxonomic_profiler == "lemur-magnet") {
            // download lemur db
            LEMUR_DATABASEDOWNLOAD(lemur_db_zenodo_id)
            LEMUR_STAGE_DB(LEMUR_DATABASEDOWNLOAD.out.db_files, LEMUR_DATABASEDOWNLOAD.out.refseq_version_bacteria)
            ch_lemur_db = LEMUR_STAGE_DB.out.lemur_db
            ch_taxonomy_db = ch_lemur_db
        } else if (params.taxonomic_profiler == "singlem") {
            SINGLEM_TAXONOMY_DOWNLOAD_DB(params.singlem_metapackage)
            ch_singlem_db = SINGLEM_TAXONOMY_DOWNLOAD_DB.out.singlem_db
            ch_taxonomy_db = ch_singlem_db
        } else if (params.taxonomic_profiler == "kraken2") {
            if (params.kraken2_db) {
                ch_kraken2_db = channel.value(file(params.kraken2_db, checkIfExists: true))
            } else {
                KRAKEN2_TAXONOMY_DOWNLOAD_DB()
                ch_kraken2_db = KRAKEN2_TAXONOMY_DOWNLOAD_DB.out.db
            }
            ch_taxonomy_db = ch_kraken2_db
        }
    }

    // ------------------------------------------------
    // assembly databases 
    // ------------------------------------------------
    if (analysis_type == "assembly" || analysis_type == "isolate-analysis") {
        // download checkm2 database 
        if (analysis_type == "assembly" || analysis_type == "isolate-analysis") {
            CHECKM2_DATABASEDOWNLOAD(checkm2_db_zenodo_id)
            ch_checkm2_db = CHECKM2_DATABASEDOWNLOAD.out.database
        }

        // download Kraken2 standard-8 database for isolate read classification
        if (analysis_type == "isolate-analysis") {
            KRAKEN2_STANDARD8_DOWNLOAD_DB()
            ch_kraken2_db = KRAKEN2_STANDARD8_DOWNLOAD_DB.out.db
        }
    
        // download bakta db
        BAKTA_BAKTADBDOWNLOAD()
        log.warn "If downloading Bakta database which is ~55GB size: it takes ~50 minutes"
        ch_bakta_db = BAKTA_BAKTADBDOWNLOAD.out.db

        // download singlem db
        if (analysis_type == "assembly") {
            SINGLEM_DOWNLOAD_DB(params.singlem_metapackage)
            ch_singlem_db = SINGLEM_DOWNLOAD_DB.out.singlem_db
        }
    }

    emit: // emit empty channels if not downloaded
    ch_deacon_index = ch_deacon_index

    // taxonomic profiling databases
    ch_emu_db = ch_emu_db // not used currently ; 
    ch_lemur_db = ch_lemur_db
    ch_sylph_db = ch_sylph_db
    
    // assembly databases
    ch_checkm2_db = ch_checkm2_db
    ch_bakta_db = ch_bakta_db
    ch_kraken2_db = ch_kraken2_db
    ch_singlem_db = ch_singlem_db
    ch_taxonomy_db = ch_taxonomy_db
}
