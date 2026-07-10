#!/usr/bin/env nextflow

include { DEACON_FETCH } from '../modules/local/deacon/fetch/main.nf'
include { DEACON_FILTER } from '../modules/local/deacon/filter/main.nf'
include { convert_to_nfcore_tuple } from '../subworkflows/local/utils/nf-core-compatibility.nf'

params.reads = "${projectDir}/../assets/data/data/metagenome_small/mock9_sub10k.fastq.gz"
params.deacon_index = 'panhuman-1'

workflow {
    reads = convert_to_nfcore_tuple(params.reads)
    DEACON_FETCH(params.deacon_index)
    DEACON_FILTER(reads, DEACON_FETCH.out.index)
}
