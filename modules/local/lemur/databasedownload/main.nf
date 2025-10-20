// note: this is based on the nf-core/checkm2 databasedownload process
def downloadZenodoApiEntry(zenodo_id) {
    // Download metadata from Zenodo API, setting "Accept: application/json" header
    def api_url  = "https://zenodo.org/api/records/${zenodo_id}"
    def conn     = new URL(api_url).openConnection()
    conn.setRequestProperty('Accept', 'application/json')
    conn.setRequestProperty('User-Agent', "Nextflow ${nextflow.version ?: ''}".trim())

    def api_text = conn.getInputStream().getText('UTF-8')
    def parser   = new groovy.json.JsonSlurper()

    return parser.parseText(api_text)
}

process LEMUR_DATABASEDOWNLOAD {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/95/95c0d3d867f5bc805b926b08ee761a993b24062739743eb82cc56363e0f7817d/data':
        'community.wave.seqera.io/library/aria2:1.37.0--3a9ec328469995dd' }"

    input:
    val(db_zenodo_id)

    output:
    val(meta), emit: db_version
    tuple path("gene2len.tsv"), 
      path("reference2genome.tsv"), 
      path("species_taxid.fasta"), 
      path("taxonomy.tsv"), emit: db_files
    val(refseq_version_bacteria) , emit: refseq_version_bacteria
    // path("versions.yml")                                 , emit: versions // not emitted since I don't have write access to the existing db directory

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    zenodo_id  = db_zenodo_id ?: 10802546  // Default to version 1 if no ID provided
    api_data   = downloadZenodoApiEntry(zenodo_id)
    db_version = api_data.metadata.version
    checksum   = api_data.files[0].checksum.replaceFirst(/^md5:/, "md5=")
    meta       = [id: 'lemur_db', version: db_version]
    db_file_name = api_data.files[0].filename ?: "rv221bacarc-rv222fungi.tar.gz" // generalize the filename for future versions and multiple files
    refseq_version_bacteria = db_file_name.find(/\d{3}/)
    """
    # download database from zenodo using aria2c (fast downloader)
    aria2c \
        ${args} \
        --checksum ${checksum} \
        https://zenodo.org/records/${zenodo_id}/files/${db_file_name}
 
    tar -xzf ${db_file_name}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        aria2: \$(echo \$(aria2c --version 2>&1) | grep 'aria2 version' | cut -f3 -d ' ')
    END_VERSIONS
    """

    stub:
    db_version = 0
    meta       = [id: 'lemur_db', version: db_version]
    """
    touch gene2len.tsv
    touch reference2genome.tsv
    touch species_taxid.fasta
    touch taxonomy.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        aria2: \$(echo \$(aria2c --version 2>&1) | grep 'aria2 version' | cut -f3 -d ' ')
    END_VERSIONS
    """
}

process LEMUR_STAGE_DB {
    label 'process_single'
    
    input:
    tuple path("gene2len.tsv"), 
      path("reference2genome.tsv"), 
      path("species_taxid.fasta"), 
      path("taxonomy.tsv")
    val(refseq_version_bacteria)
    
    output:
    path "lemur_${refseq_version_bacteria}_db/", emit: lemur_db
    
    script:
    """
    mkdir lemur_${refseq_version_bacteria}_db

    mv *.tsv lemur_${refseq_version_bacteria}_db/
    mv *.fasta lemur_${refseq_version_bacteria}_db/
    """
    
    stub:
    """
    mkdir lemur_221_db
    mv *.tsv lemur_221_db/
    mv *.fasta lemur_221_db/
    """
}
