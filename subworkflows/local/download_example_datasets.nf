#!/usr/bin/env nextflow

params.example_datasets_url           = "https://drive.google.com/drive/u/1/folders/1vPR4hEGiYfeoG4Vy1qQ6RuD0UqlYhQAH"
// params.example_datasets_url = "https://drive.google.com/drive/folders/1HRYK6EW9HVThJfeNObvMo_fCo0Dp_gkA?usp=sharing" // dummy url for testing


process GDOWN {
    label 'process_single'
    conda "conda-forge::gdown"

    input:
    val google_drive_url
    val save_dir

    output:
    path "download_complete.txt"

    script:
    """
    mkdir -p "${save_dir}"
    gdown --folder "${google_drive_url}" -c -O "${save_dir}"
    find "${save_dir}" -type f | sort > download_complete.txt
    """
}


workflow {
    save_dir = params.containsKey('save_dir') && params.save_dir
        ? params.save_dir
        : (params.containsKey('outdir') && params.outdir ? params.outdir : "${launchDir}/assets/data")

    GDOWN(params.example_datasets_url, save_dir)
}
