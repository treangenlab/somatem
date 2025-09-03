#!/usr/bin/env nextflow

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

workflow TAXONOMIC_PROFILING {

    take:
    clean_reads_ch

   main:
    // 16S amplicon reads
    if (params.data_type == "16S") {
        EMU_ABUNDANCE(clean_reads_ch)  

    // metagenomic reads (default)
    } else {
        LEMUR(clean_reads_ch) // tax profiling

        lemur_classification = LEMUR.out.output_dir.map { dir -> dir + '/relative_abundance.tsv' } // gather classification file
        MAGNET(clean_reads_ch, lemur_classification) // Correct false positives

    }

    // emit:
    // abundance = MAGNET.out // TODO: use if/else to emit different channels?
}
