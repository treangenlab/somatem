#!/usr/bin/env nextflow

include { convert_to_nfcore_tuple } from './utils/nf-core-compatibility.nf'

include { EMU_ABUNDANCE } from '../../modules/local/emu/main.nf'
include { LEMUR } from '../../modules/local/lemur/main.nf'
include { MAGNET } from '../../modules/local/magnet/main.nf'




// -------------------------
// Parameters
// -------------------------

params.input_dir   = 'examples/data'

// other params
params.lemur_db = "${projectDir}/../examples/lemur/example-db/"
params.lemur_taxonomy = "${projectDir}/../examples/lemur/example-db/taxonomy.tsv"
params.rank = "species"

// -------------------------
// Workflow Definition
// -------------------------

workflow {
    reads_ch = convert_to_nfcore_tuple(params.input_dir)
    
    // 16S amplicon reads
    if (params.data_type == "16S") {
        EMU_ABUNDANCE(reads_ch)  

    // metagenomic reads (default)
    } else {
        database_dir = Channel.fromPath(params.lemur_db)
        taxonomy = Channel.fromPath(params.lemur_taxonomy)
        rank = Channel.of(params.rank)

        LEMUR(reads_ch, database_dir, taxonomy, rank) // tax profiling

        lemur_classification = LEMUR.out.map { dir -> dir + '/relative_abundance.tsv' }
        MAGNET(reads_ch, lemur_classification) // Correct false positives

    }
}
