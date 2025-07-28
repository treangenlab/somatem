// script for testing nextflow syntax and stuff

params.input_dir = 'examples/data/'

include { convert_to_nfcore_tuple } from './subworkflows/utils/nf-core-compatibility.nf'

workflow {

// in_ch = Channel.fromPath("${params.input_dir}/*.fastq.gz")
// in_ch.view { r -> "files: ${r.simpleName}"}

// test the convert_to_nfcore_tuple subworkflow
in_ch = convert_to_nfcore_tuple(params.input_dir)
in_ch.view { r -> "tuple: ${r}"}

// // test using queue channel twice
// print_name(in_ch) | view()

// // in_ch.view { r -> "used: ${r}"}

// print_extension(in_ch) | view()

// // test empty channels
// contam_ref = Channel.of([])
// temp_ch = Channel.empty()

// contam_ref.view { r -> "empty: ${r}"}
// temp_ch.view { r -> "empty: ${r}"}

}


process print_name {
    input:
    tuple val(meta), path(reads)
    output:
    stdout    
    script:
    """
    echo "${reads.simpleName}"
    """
}

process print_extension {
    input:
    tuple val(meta), path(reads)
    output:
    stdout    
    script:
    """
    echo "${reads.extension}"
    """
}