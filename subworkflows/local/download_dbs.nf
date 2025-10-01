#!/usr/bin/env nextflow

// include { checkm2_DOWNLOAD_DB } from '../../modules/local/checkm2/download_db/main.nf'
include { CHECKM2_DATABASEDOWNLOAD } from '../../modules/nf-core/checkm2/databasedownload/main'   
include { BAKTA_BAKTADBDOWNLOAD } from '../../modules/nf-core/bakta/baktadbdownload/main' 
include { SINGLEM_DOWNLOAD_DB } from '../../modules/local/singlem/download_db/main.nf'

workflow DOWNLOAD_DBS {

    take:
    None

    main:
    ch_dbs = Channel.empty()

    // download checkm2 database 
    CHECKM2_DATABASEDOWNLOAD(params.checkm2_db_zenodo_id)
    ch_dbs = ch_dbs.mix(CHECKM2_DATABASEDOWNLOAD.out.database)
    
    // download bakta db
    BAKTA_BAKTADBDOWNLOAD()
    log.warn "If downloading Bakta database which is ~55GB size: it takes ~50 minutes"
    ch_dbs = ch_dbs.mix(BAKTA_BAKTADBDOWNLOAD.out.db)

    // download singlem db
    SINGLEM_DOWNLOAD_DB()
    ch_dbs = ch_dbs.mix(SINGLEM_DOWNLOAD_DB.out.singlem_db)

    emit:
    dbs = ch_dbs
}
