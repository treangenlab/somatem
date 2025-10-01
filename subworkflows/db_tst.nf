#!/usr/bin/env nextflow

include { DOWNLOAD_DBS } from './local/download_dbs'

workflow {

    ch_dbs = Channel.empty()

    // download checkm2 database 
    DOWNLOAD_DBS(0)
    ch_dbs = ch_dbs.mix(DOWNLOAD_DBS.out.dbs)
    
    // print dbs
    ch_dbs.view { x -> "Input reads: ${x}" }
}