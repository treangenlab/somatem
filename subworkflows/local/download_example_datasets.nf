#!/usr/bin/env nextflow
// include { GDOWN } from "../../modules/local/gdown/main"

process GDOWN2 {
    label 'process_single'
    conda "conda-forge::gdown"

    input:
    val google_drive_url
    path save_dir

    script:
    """
    # gdown --folder ${google_drive_url} -c -O ${save_dir}
    echo "test parameter: data_type = ${params.data_type}"
    """
}


workflow DOWNLOAD_EXAMPLE_DATASETS {

    take:
    google_drive_url
    save_dir

    main:
    GDOWN2(google_drive_url, save_dir)
}