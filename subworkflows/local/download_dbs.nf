#!/usr/bin/env nextflow

// include { checkm2_DOWNLOAD_DB } from '../../modules/local/checkm2/download_db/main.nf'
include { CHECKM2_DATABASEDOWNLOAD } from '../../modules/nf-core/checkm2/databasedownload/main'   
include { BAKTA_BAKTADBDOWNLOAD } from '../../modules/nf-core/bakta/baktadbdownload/main' 
include { SINGLEM_DOWNLOAD_DB } from '../../modules/local/singlem/download_db/main.nf'

workflow {

    // download checkm2 database 
    CHECKM2_DATABASEDOWNLOAD(params.checkm2_db_zenodo_id)
    log.info "Downloading CheckM2 database, 3GB size: Takes ~3 minutes"

    
    // download bakta db
    BAKTA_BAKTADBDOWNLOAD()
    log.info "Downloading Bakta database, ~36GB size: Takes >30 minutes"

    // download singlem db
    SINGLEM_DOWNLOAD_DB()
    log.info "Downloading SingleM database, <2 GB size: Takes ~3 minutes"

}
