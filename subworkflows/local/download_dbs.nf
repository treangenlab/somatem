#!/usr/bin/env nextflow

// include { checkm2_DOWNLOAD_DB } from '../../modules/local/checkm2/download_db/main.nf'
include { CHECKM2_DATABASEDOWNLOAD } from '../../modules/nf-core/checkm2/databasedownload/main'   
include { BAKTA_BAKTADBDOWNLOAD } from '../../modules/nf-core/bakta/baktadbdownload/main' 
include { SINGLEM_DOWNLOAD_DB } from '../../modules/local/singlem/download_db/main.nf'

workflow {

    // download checkm2 database 
    CHECKM2_DATABASEDOWNLOAD(params.checkm2_db_zenodo_id)
    
    // download bakta db
    BAKTA_BAKTADBDOWNLOAD()
    log.warn "If downloading Bakta database which is ~55GB size: it takes ~50 minutes"

    // download singlem db
    SINGLEM_DOWNLOAD_DB()

}
