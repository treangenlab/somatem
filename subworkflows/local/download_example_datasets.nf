#!/usr/bin/env nextflow

params.example_datasets_zenodo_record = '21299097'
params.example_datasets_url = 'https://drive.google.com/drive/u/1/folders/1vPR4hEGiYfeoG4Vy1qQ6RuD0UqlYhQAH'

process DOWNLOAD_ZENODO_EXAMPLE_DATASETS {
    label 'process_single'
    conda 'conda-forge::python>=3.11 conda-forge::gdown'

    input:
    val record_id
    val google_drive_url
    val save_dir

    output:
    path 'download_complete.txt'

    script:
    """
    if ! python ${projectDir}/../../bin/download_zenodo_examples.py \\
            --record-id '${record_id}' \\
            --save-dir '${save_dir}' \\
            --manifest download_complete.txt; then
        echo 'Zenodo download failed; trying the Google Drive backup.' >&2
        mkdir -p '${save_dir}'
        gdown --folder '${google_drive_url}' -c -O '${save_dir}'
        python ${projectDir}/../../bin/download_zenodo_examples.py \
            --record-id '${record_id}' \
            --save-dir '${save_dir}' \
            --manifest download_complete.txt \
            --samplesheets-only
        find '${save_dir}' -type f | sort > download_complete.txt
    fi
    """
}

workflow {
    save_dir = params.containsKey('save_dir') && params.save_dir
        ? params.save_dir
        : (params.containsKey('outdir') && params.outdir ? params.outdir : "${launchDir}/assets/data")
    save_dir = file(save_dir).toAbsolutePath().toString()

    DOWNLOAD_ZENODO_EXAMPLE_DATASETS(
        params.example_datasets_zenodo_record,
        params.example_datasets_url,
        save_dir
    )
}
