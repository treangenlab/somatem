#!/usr/bin/env nextflow

// entry workflow for downloading databases before starting the main workflow
include { DOWNLOAD_DBS } from '../download_databases.nf'

workflow {

    DOWNLOAD_DBS(params.analysis_type, params.hostile_index, 
            params.lemur_db_zenodo_id, params.checkm2_db_zenodo_id, params.target_taxonomic_domain, params.sylph_db_picker)
}