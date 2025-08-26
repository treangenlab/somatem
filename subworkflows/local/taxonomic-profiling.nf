#!/usr/bin/env nextflow

include { convert_to_nfcore_tuple } from './utils/nf-core-compatibility.nf'

include { EMU_ABUNDANCE } from '../../modules/local/emu/main.nf'
include { LEMUR } from '../../modules/local/lemur/main.nf'
include { MAGNET } from '../../modules/local/magnet/main.nf'

// -------------------------
// Parameters
// -------------------------
params.input_dir   = 'examples/data'


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
        LEMUR(reads_ch) // tax profiling

        lemur_classification = LEMUR.out.output_dir.map { dir -> dir + '/relative_abundance.tsv' } // gather classification file
        MAGNET(reads_ch, lemur_classification) // Correct false positives

    }
}
