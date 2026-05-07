#!/usr/bin/env nextflow

params.example_datasets_url           = "https://drive.google.com/drive/u/1/folders/1vPR4hEGiYfeoG4Vy1qQ6RuD0UqlYhQAH"
// params.example_datasets_url = "https://drive.google.com/drive/folders/1HRYK6EW9HVThJfeNObvMo_fCo0Dp_gkA?usp=sharing" // dummy url for testing

params.save_dir = params.save_dir ?: "${launchDir}/assets" // default save directory: 
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