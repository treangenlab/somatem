#!/usr/bin/env nextflow

params.example_datasets_url           = "https://drive.google.com/drive/folders/11ZRpUCRrhdcJarlYdMSEDlCFl3oIz6Bh?usp=sharing"
params.save_dir = params.save_dir ?: "${launchDir}/somatem_assets" // default save directory: 
  // currently users can't change this since the config yml files are hardcoded to use this path


process GDOWN {
    label 'process_single'
    conda "conda-forge::gdown"

    input:
    val google_drive_url
    val save_dir

    script:
    """
    gdown --folder ${google_drive_url} -c -O ${save_dir}
    """
}


workflow {
    GDOWN(params.example_datasets_url, params.save_dir)
}