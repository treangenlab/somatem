process KRAKEN2_STANDARD8_DOWNLOAD_DB {
    tag 'kraken2_standard8_download_db'
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/95/95c0d3d867f5bc805b926b08ee761a993b24062739743eb82cc56363e0f7817d/data' :
        'community.wave.seqera.io/library/aria2:1.37.0--3a9ec328469995dd' }"

    output:
    path "k2_standard_08_GB_20260226", emit: db
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def db_url = 'https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08_GB_20260226.tar.gz'
    def db_archive = db_url.tokenize('/')[-1]
    def db_dir = db_archive.replaceFirst(/\.tar\.gz$/, '')
    """
    aria2c \\
        $args \\
        --out ${db_archive} \\
        ${db_url}

    mkdir -p ${db_dir} db_tmp
    tar -xzf ${db_archive} -C db_tmp

    if [ -f db_tmp/hash.k2d ]; then
        mv db_tmp/* ${db_dir}/
    else
        extracted_db_dir=\$(find db_tmp -name hash.k2d -printf '%h\\n' | head -n 1)
        if [ -z "\${extracted_db_dir}" ]; then
            echo "Could not find hash.k2d in extracted Kraken2 database archive" >&2
            exit 1
        fi
        mv "\${extracted_db_dir}"/* ${db_dir}/
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        aria2: \$(echo \$(aria2c --version 2>&1) | grep 'aria2 version' | cut -f3 -d ' ')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p k2_standard_08_GB_20260226
    touch k2_standard_08_GB_20260226/hash.k2d
    touch k2_standard_08_GB_20260226/opts.k2d
    touch k2_standard_08_GB_20260226/taxo.k2d

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        aria2: "stub"
    END_VERSIONS
    """
}
