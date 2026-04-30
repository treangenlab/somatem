// downloads test data from google drive `examples` directory 

process GDOWN {
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    input:
    val google_drive_url
    path save_dir

    script:
    """
    echo "Downloading from: ${google_drive_url}"
    gdown --folder ${google_drive_url} -c -O ${save_dir}
    """
}

