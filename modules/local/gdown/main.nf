// downloads test data from google drive `examples` directory 

process GDOWN {
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    script:
    base_dir = "${moduleDir}/../../.."
    """
    gdown --folder https://drive.google.com/drive/folders/11ZRpUCRrhdcJarlYdMSEDlCFl3oIz6Bh?usp=sharing -c -O ${base_dir}/assets/
    """
}

