#!/usr/bin/env nextflow

include { checkm2_DOWNLOAD_DB } from '../../modules/local/checkm2/download_db/main.nf'
    

workflow {
    checkm2_DOWNLOAD_DB(params.checkm2_db)
}
