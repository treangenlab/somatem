#!/usr/bin/env nextflow

// include { GDOWN } from '../../modules/local/gdown/main.nf'

params.outdir = params.outdir ?: "${launchDir}/somatem_example_data"
params.example_datasets_url = params.example_datasets_url ?: "https://drive.google.com/drive/folders/1HRYK6EW9HVThJfeNObvMo_fCo0Dp_gkA?usp=sharing" // dummy url for testing
params.save_dir = params.save_dir ?: "${launchDir}/somatem_assets" // default save directory


process GDOWN2 {
    label 'process_single'
    conda "conda-forge::gdown"

    input:
    val google_drive_url
    path save_dir

    script:
    """
    gdown --folder ${google_drive_url} -c -O ${save_dir}
    echo "test parameter: save_dir = ${params.save_dir}"
    echo "launchDir = ${launchDir}"
    """
}


workflow {
    GDOWN2(params.example_datasets_url, params.save_dir)
}