#!/usr/bin/env nextflow

// include { checkm2_DOWNLOAD_DB } from '../../modules/local/checkm2/download_db/main.nf'
include { CHECKM2_DATABASEDOWNLOAD } from '../../modules/nf-core/checkm2/databasedownload/main'    

workflow {

    // download checkm2 database 
    CHECKM2_DATABASEDOWNLOAD(params.checkm2_db_zenodo_id)
    
}
