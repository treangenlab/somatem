/*
 * MARLIN COMPLETE SUBWORKFLOW
 * 
 * Complete implementation of the Marlin pipeline including:
 * - Lemur taxonomic classification
 * - MAGnet reference selection based on taxonomy
 * - Single-pass assembly with adaptive variant calling (using MAGnet-selected refs)
 * - Contig breaking at coverage gaps
 * - dRep genome dereplication
 * 
 * Workflow order:
 * 1. Lemur classifies reads taxonomically
 * 2. MAGnet uses classifications to select relevant references
 * 3. For each read/reference pair:
 *    a. Align reads with Minimap2
 *    b. Calculate depth coverage
 *    c. Call variants adaptively (Clair3 or bcftools based on coverage)
 *    d. Classify variants by confidence (HIGH/MED/LOW)
 *    e. Generate full-length consensus from HIGH-confidence variants
 *    f. Break consensus at coverage gaps to create contigs
 * 4. dRep dereplicates the resulting contigs
 */

// Marlin single-pass assembly modules
include { MINIMAP2_ALIGN     } from '../../modules/nf-core/minimap2/align'
include { SAMTOOLS_DEPTH     } from '../../modules/nf-core/samtools/depth'
include { CLAIR3             } from '../../modules/local/clair3'
include { CLASSIFY_VARIANTS  } from '../../modules/local/marlin/classify_variants'
include { BCFTOOLS_CONSENSUS } from '../../modules/local/bcftools/consensus'
include { BUILD_CONTIGS      } from '../../modules/local/marlin/build_contigs'

// Taxonomic and reference selection modules
include { LEMUR                   } from '../../modules/local/lemur'
include { SUGGEST_MAGNET_CUTOFF   } from '../../modules/local/marlin/magnet_cutoff'
include { MAGNET                  } from '../../modules/local/magnet'
include { PARSE_MAGNET_REFS       } from '../../modules/local/marlin/parse_refs'

// Adaptive variant calling modules
include { SELECT_VARIANT_CALLER   } from '../../modules/local/marlin/select_caller'
include { BCFTOOLS_CALL           } from '../../modules/local/bcftools/call'

// Dereplication module
include { DREP_DEREPLICATE        } from '../../modules/local/drep'

workflow MARLIN_COMPLETE {
    take:
    ch_reads            // channel: [ val(meta), path(reads) ]
    ch_lemur_db         // channel: path(lemur_database)
    clair3_model        // val: path to Clair3 model (optional)
    clair3_platform     // val: Clair3 platform (default: ont)
    ani_threshold       // val: ANI threshold for dRep (default: 0.95)
    af_high             // val: AF threshold for HIGH confidence (default: 0.9)
    af_med              // val: AF threshold for MED confidence (default: 0.7)
    high_cov_frac       // val: Fraction of median for high coverage (default: 0.5)
    gap_min_cov_frac    // val: Fraction of median for gap detection (default: 0.1)
    gap_min_cov_floor   // val: Minimum coverage floor (default: 2)
    gap_min_len         // val: Minimum gap length to break contigs (default: 200)
    contig_min_len      // val: Minimum contig length to output (default: 1000)
    magnet_target_frac  // val: Target fraction for MAGnet cutoff (default: 0.99)
    min_cov_for_clair3  // val: Minimum coverage to use Clair3 vs bcftools (default: 8)

    main:
    ch_versions = channel.empty()

    //
    // STEP 1: Run Lemur taxonomic classification on input reads
    //
    LEMUR(
        ch_reads,
        ch_lemur_db
    )
    ch_versions = ch_versions.mix(LEMUR.out.versions)

    //
    // STEP 1b: Suggest optimal MAGnet cutoff from Lemur abundance table
    //
    // Analyzes Lemur abundance profile to automatically determine --min-abundance
    // cutoff that captures the target fraction (default 99%) of community abundance
    SUGGEST_MAGNET_CUTOFF(
        LEMUR.out.report,
        magnet_target_frac
    )
    ch_versions = ch_versions.mix(SUGGEST_MAGNET_CUTOFF.out.versions)

    //
    // STEP 2: Run MAGnet to select relevant references based on Lemur classifications
    //
    // MAGnet takes reads + Lemur classifications and outputs which references to use
    // The auto-calculated cutoff is used via task.ext.args in nextflow.config
    ch_magnet_input = ch_reads
        .join(LEMUR.out.report)

    MAGNET(
        ch_magnet_input.map { meta, reads, _classification -> [meta, reads] },
        ch_magnet_input.map { _meta, _reads, classification -> classification }
    )
    ch_versions = ch_versions.mix(MAGNET.out.versions)

    //
    // STEP 3: Parse MAGnet results to extract "Present" references
    //
    // MAGnet.out.report (cluster_representative.csv) contains presence/absence calls
    // Extract only references marked as "Present" or "Genus Present"
    PARSE_MAGNET_REFS(
        MAGNET.out.report
    )
    ch_versions = ch_versions.mix(PARSE_MAGNET_REFS.out.versions)

    //
    // STEP 4: Extract reference FASTAs from MAGnet output directory
    //
    // MAGnet outputs a reference_genomes/ directory containing FASTA files
    ch_magnet_refs = MAGNET.out.ref_dir
        .map { meta, ref_dir ->
            // Get all FASTA files from reference_genomes directory
            def ref_files = []
            def dir_file = ref_dir instanceof List ? ref_dir[0] : ref_dir
            
            if (dir_file.isDirectory()) {
                dir_file.eachFileMatch(~/.*\.(fa|fasta|fna)$/) { file ->
                    ref_files << file
                }
            }
            
            [ meta, ref_files ]
        }
        .transpose()  // Flatten to one reference per channel item

    //
    // STEP 5: Filter MAGnet references to only "Present" ones
    //
    // Join the parsed reference IDs with the MAGnet-downloaded reference files
    // Only keep references that are in the "Present" list
    ch_filtered_refs = ch_magnet_refs
        .combine(PARSE_MAGNET_REFS.out.ref_ids, by: 0)
        .map { meta, ref_fasta, ref_ids_file ->
            // Extract ref ID from filename
            def ref_id = ref_fasta.getBaseName()
            [ meta, ref_id, ref_fasta, ref_ids_file ]
        }
        .filter { _meta, ref_id, _ref_fasta, ref_ids_file ->
            // Check if this reference ID is in the selected list
            def selected_ids = ref_ids_file.text.readLines()*.trim()
            selected_ids.contains(ref_id)
        }
        .map { meta, _ref_id, ref_fasta, _ref_ids_file ->
            [ meta, ref_fasta ]
        }

    //
    // STEP 6: Create read-reference pairs using filtered MAGnet references
    //
    // Combine each sample's reads with each of its MAGnet-selected "Present" references
    ch_read_ref_pairs = ch_reads
        .join(ch_filtered_refs.groupTuple())  // Group all selected refs per sample
        .transpose()                           // Create one pair per reference
        .map { meta, reads, ref_fasta ->
            // Extract reference ID from filename
            def ref_id = ref_fasta.getBaseName()
            def combined_meta = [
                id: "${meta.id}_vs_${ref_id}",
                sample_id: meta.id,
                ref_id: ref_id
            ]
            [ combined_meta, reads, ref_fasta ]
        }

    //
    // STEP 7: Run single-pass assembly for MAGnet-selected read-reference combinations
    //

    //
    // STEP 7a: Align reads to reference with Minimap2
    //
    MINIMAP2_ALIGN(
        ch_read_ref_pairs.map { meta, reads, ref_fasta -> [ meta, reads ] },
        ch_read_ref_pairs.map { meta, reads, ref_fasta -> [ meta, ref_fasta ] },
        true,      // bam_format - output as BAM
        'bai',     // bam_index_extension
        false,     // cigar_paf_format
        false      // cigar_bam
    )
    ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions)

    //
    // STEP 7b: Calculate depth coverage (needed for variant caller selection)
    //
    // SAMTOOLS_DEPTH expects: tuple val(meta), path(bam), path(bai), path(reference)
    // Join BAM + BAI + reference FASTA
    ch_depth_input = MINIMAP2_ALIGN.out.bam
        .join(MINIMAP2_ALIGN.out.index, by: 0)
        .join(ch_read_ref_pairs.map { meta, reads, ref_fasta -> [ meta, ref_fasta ] }, by: 0)
    
    SAMTOOLS_DEPTH(
        ch_depth_input
    )
    ch_versions = ch_versions.mix(SAMTOOLS_DEPTH.out.versions)

    //
    // STEP 7c: Select variant caller based on coverage
    //
    // Analyze coverage to determine whether to use Clair3 (≥8x) or bcftools (<8x)
    SELECT_VARIANT_CALLER(
        SAMTOOLS_DEPTH.out.depth,
        min_cov_for_clair3
    )
    ch_versions = ch_versions.mix(SELECT_VARIANT_CALLER.out.versions)

    //
    // STEP 7d: Branch based on variant caller choice
    //
    ch_caller_choice = SELECT_VARIANT_CALLER.out.choice
        .map { meta, choice_file ->
            def caller = choice_file.text.trim()
            [ meta, caller ]
        }

    // Split BAM channel based on caller choice
    // MINIMAP2_ALIGN.out.bam emits [meta, bam] (2 elements)
    // We need to reconstruct the full tuple with BAI and reference
    ch_bam_with_ref = MINIMAP2_ALIGN.out.bam
        .join(ch_read_ref_pairs.map { meta, reads, ref_fasta -> [ meta, ref_fasta ] })
        .map { meta, bam, ref_fasta ->
            def bai = file("${bam}.bai")
            [ meta, bam, bai, ref_fasta ]
        }
    
    ch_bam_with_choice = ch_bam_with_ref
        .join(ch_caller_choice)
        .branch { meta, bam, bai, ref, caller ->
            clair3: caller == 'clair3'
                return [ meta, bam, bai, ref ]
            bcftools: caller == 'bcftools'
                return [ meta, bam, bai, ref ]
        }

    //
    // STEP 7e: Call variants with Clair3 (high coverage samples)
    //
    CLAIR3(
        ch_bam_with_choice.clair3,
        [],
        clair3_platform
    )
    ch_versions = ch_versions.mix(CLAIR3.out.versions)

    //
    // STEP 7f: Call variants with bcftools (low coverage samples)
    //
    BCFTOOLS_CALL(
        ch_bam_with_choice.bcftools
    )
    ch_versions = ch_versions.mix(BCFTOOLS_CALL.out.versions)

    //
    // STEP 7g: Merge VCF outputs from both variant callers
    //
    ch_all_vcfs = CLAIR3.out.vcf.mix(BCFTOOLS_CALL.out.vcf)

    //
    // STEP 7h: Classify variants using depth and allele frequency
    //
    // Join VCF with depth file for classification
    ch_classify_input = ch_all_vcfs
        .join(SAMTOOLS_DEPTH.out.depth)
        .map { meta, vcf, tbi, depth_file ->
            [ meta, vcf, tbi, depth_file ]
        }

    CLASSIFY_VARIANTS(
        ch_classify_input,
        af_high,
        af_med,
        high_cov_frac,
        gap_min_cov_frac,
        gap_min_cov_floor
    )
    ch_versions = ch_versions.mix(CLASSIFY_VARIANTS.out.versions)

    //
    // STEP 7i: Generate consensus sequence with BCFtools using HIGH-confidence variants
    //
    // Join reference with high-confidence VCF for consensus generation
    // Get reference from ch_read_ref_pairs since MINIMAP2_ALIGN.out.bam doesn't include it
    ch_reference_map = ch_read_ref_pairs.map { meta, reads, ref_fasta -> [ meta, ref_fasta ] }
    
    ch_consensus_input = ch_reference_map
        .join(CLASSIFY_VARIANTS.out.highconf_vcf)
        .map { meta, ref, vcf, tbi ->
            [ meta, ref, vcf, tbi ]
        }

    BCFTOOLS_CONSENSUS(
        ch_consensus_input
    )
    ch_versions = ch_versions.mix(BCFTOOLS_CONSENSUS.out.versions)

    //
    // STEP 7j: Break consensus at coverage gaps to create contigs
    //
    // Join consensus FASTA with depth file for gap-based contig breaking
    // min_cov is calculated as: max(gap_min_cov_floor, global_median * gap_min_cov_frac)
    // This matches the Python script's behavior
    ch_build_contigs_input = BCFTOOLS_CONSENSUS.out.fasta
        .join(SAMTOOLS_DEPTH.out.depth)
        .map { meta, consensus, depth_file ->
            [ meta, consensus, depth_file ]
        }

    BUILD_CONTIGS(
        ch_build_contigs_input,
        gap_min_cov_floor,  // Use floor as min_cov (script will compute actual threshold)
        gap_min_len,
        contig_min_len
    )
    ch_versions = ch_versions.mix(BUILD_CONTIGS.out.versions)

    //
    // STEP 8: Run dRep dereplication on assembled genomes
    //
    // Group all contig FASTAs per sample for dereplication
    // Use BUILD_CONTIGS output (consensus_contigs.fa) instead of full consensus
    ch_drep_input = BUILD_CONTIGS.out.contigs
        .map { meta, contigs ->
            [ [id: meta.sample_id], contigs ]
        }
        .groupTuple()

    DREP_DEREPLICATE(
        ch_drep_input,
        ani_threshold
    )
    ch_versions = ch_versions.mix(DREP_DEREPLICATE.out.versions)

    emit:
    // Lemur outputs
    lemur_classifications = LEMUR.out.report                        // [ meta, tsv ]
    
    // MAGnet cutoff outputs
    magnet_cutoff         = SUGGEST_MAGNET_CUTOFF.out.cutoff        // [ meta, txt ]
    magnet_cutoff_stats   = SUGGEST_MAGNET_CUTOFF.out.stats         // [ meta, tsv ]
    
    // MAGnet outputs (reference selection and analysis)
    magnet_report         = MAGNET.out.report                       // [ meta, csv ]
    magnet_refs           = MAGNET.out.ref_dir                      // [ meta, dir ]
    magnet_bam            = MAGNET.out.bam                          // [ meta, bam ]
    magnet_selected_ids   = PARSE_MAGNET_REFS.out.ref_ids          // [ meta, txt ]
    
    // Variant caller selection outputs
    variant_caller_choice = SELECT_VARIANT_CALLER.out.choice        // [ meta, txt ]
    coverage_stats        = SELECT_VARIANT_CALLER.out.stats         // [ meta, tsv ]
    
    // Assembly outputs (only for MAGnet-selected references)
    consensus_full        = BCFTOOLS_CONSENSUS.out.fasta            // [ meta, fasta ]
    consensus_contigs     = BUILD_CONTIGS.out.contigs               // [ meta, fasta ]
    contig_layout         = BUILD_CONTIGS.out.layout                // [ meta, tsv ]
    contigs               = BUILD_CONTIGS.out.contigs               // [ meta, fasta ]
    variants              = ch_all_vcfs                             // [ meta, vcf.gz, tbi ]
    classified_variants   = CLASSIFY_VARIANTS.out.classified_vcf    // [ meta, vcf ]
    highconf_variants     = CLASSIFY_VARIANTS.out.highconf_vcf      // [ meta, vcf.gz, tbi ]
    alignments            = MINIMAP2_ALIGN.out.bam                  // [ meta, bam ]
    depth                 = SAMTOOLS_DEPTH.out.depth                // [ meta, txt ]
    
    // dRep outputs
    dereplicated_genomes  = DREP_DEREPLICATE.out.genomes            // [ meta, fasta ]
    drep_tables           = DREP_DEREPLICATE.out.tables             // [ meta, csv ]
    drep_figures          = DREP_DEREPLICATE.out.figures            // [ meta, pdf ]
    drep_logs             = DREP_DEREPLICATE.out.logs               // [ meta, log ]
    
    versions              = ch_versions
}

/*
 * Entry workflow for standalone execution
 */
workflow {
    // Validate required parameters
    if (!params.reads) {
        error "Missing required parameter: --reads"
    }
    if (!params.lemur_db) {
        error "Missing required parameter: --lemur_db"
    }

    // Create input channels
    ch_reads = channel.fromPath(params.reads)
        .map { file ->
            def meta = [id: file.getSimpleName()]
            [meta, file]
        }
    
    ch_lemur_db = channel.fromPath(params.lemur_db)

    // Run the subworkflow
    MARLIN_COMPLETE(
        ch_reads,
        ch_lemur_db,
        params.clair3_model ?: null,
        params.clair3_platform ?: 'ont',
        params.ani_threshold ?: 0.95,
        params.af_high ?: 0.9,
        params.af_med ?: 0.7,
        params.high_cov_frac ?: 0.5,
        params.gap_min_cov_frac ?: 0.1,
        params.gap_min_cov_floor ?: 2,
        params.gap_min_len ?: 200,
        params.contig_min_len ?: 1000,
        params.magnet_target_frac ?: 0.99,
        params.min_cov_for_clair3 ?: 8
    )

    // View key outputs
    MARLIN_COMPLETE.out.dereplicated_genomes.view { meta, fasta ->
        "Sample ${meta.id}: Dereplicated genomes -> ${fasta}"
    }
}
