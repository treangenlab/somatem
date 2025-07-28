#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// -------------------------
// Parameters (all configurable)
// -------------------------
params.input_dir = '/path/to/your/fastq/folder'
params.out_dir = '/home/tmhagm8/Documents/SOMAteM/examples/soma_test_out'
params.threads = 16
params.gtdb_db = '/path/to/gtdb/database'

// Note: out_dir normalization moved to workflow to avoid redefinition warning


// -------------------------
// Process: myloasm (assembly & polish)
//    - Takes (sample_id, reads) → emits (sample_id, assembly_primary.fa)
// -------------------------
process runMyloasm {
    tag "${sample_id}"
    cpus params.threads
    conda 'bioconda::myloasm'

    publishDir path: "${params.out_dir}/myloasm", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("assembly_primary.fa"), emit: polished_tuple

    script:
    // Join reads into one string if it's a list
    def reads_list = (reads instanceof List) ? reads.join(' ') : reads.toString()
    """
    myloasm ${reads_list} \
      -o . \
      -t ${task.cpus}
    """
}


// -------------------------
// Process: runSingleM (taxonomic profiling)
//    - Takes (sample_id, polished_fasta) → emits (sample_id, singlem_report.tsv)
// -------------------------
process runSingleM {
    tag "${sample_id}"
    cpus params.threads
    conda 'bioconda::singlem'

    publishDir path: "${params.out_dir}/singlem", mode: 'copy'

    input:
    tuple val(sample_id), path(polished_fasta)

    output:
    tuple val(sample_id), path("singlem_report.tsv"), emit: singlem_tuple

    script:
    """
    singlem pipe \\
      --genome-fasta-files ${polished_fasta} \\
      --taxonomic-profile singlem_report.tsv \\
      --threads ${task.cpus}
    """
}


// -------------------------
// Process: runCoverM (coverage of contigs vs. reads)
//    - Takes (sample_id, polished_fasta, reads) → emits (sample_id, depth.tsv)
// -------------------------
process runCoverM {
    tag "${sample_id}"
    cpus params.threads
    conda 'bioconda::coverm'

    publishDir path: "${params.out_dir}/coverm", mode: 'copy'

    input:
    tuple val(sample_id), path(polished_fasta), path(reads)

    output:
    tuple val(sample_id), path("depth.tsv"), emit: depth_tuple

    script:
    def reads_list = (reads instanceof List) ? reads.join(' ') : reads.toString()
    """
    coverm contig \
      --coupled ${reads_list} \
      --methods metabat \
      --reference ${polished_fasta} \
      --threads ${task.cpus} \
      --output-file depth.tsv
    """
}


// -------------------------
// Process: runSemibin2 (binning #1)
//    - Takes (sample_id, polished_fasta, depth_file) → emits (sample_id, semibin2_bins/)
// -------------------------
process runSemibin2 {
    tag "${sample_id}"
    cpus params.threads
    conda 'bioconda::semibin'

    publishDir path: "${params.out_dir}/semibin2", mode: 'copy'

    input:
    tuple val(sample_id), path(polished_fasta), path(depth_file)

    output:
    tuple val(sample_id), path("semibin2_bins/semibin2_bin.*.fna"), emit: semibin_bins

    script:
    """
    mkdir -p semibin2_bins
    semibin2 single_easy_bin \
      --input-fasta ${polished_fasta} \
      --depth-file ${depth_file} \
      --output semibin2_bins \
      --threads ${task.cpus}
    """
}


// -------------------------
// Process: runMetabat2 (binning #2)
//    - Takes (sample_id, polished_fasta, depth_file) → emits (sample_id, metabat2_bins/)
// -------------------------
process runMetabat2 {
    tag { sample_id }
    cpus params.threads
    conda 'bioconda::metabat2'

    publishDir path: "${params.out_dir}/metabat2", mode: 'copy'

    input:
    tuple val(sample_id), path(polished_fasta), path(depth_file)

    output:
    tuple val(sample_id), path("metabat2_bins/metabat2_bin.*.fna"), emit: metabat2_bins

    script:
    """
    # Create output directory for bins
    mkdir -p metabat2_bins

    # Run MetaBAT2 with precomputed depth file
    metabat2 \
      -i ${polished_fasta} \
      -a ${depth_file} \
      -t ${task.cpus} \
      -m 1500 \
      -o metabat2_bins/bin \
      --verbose
    """
}


// -------------------------
// Process: runVamb (binning #3)
//    - Takes (sample_id, polished_fasta, depth_file) → emits (sample_id, vamb_bins/)
// -------------------------
process runVamb {
    tag { sample_id }
    cpus params.threads
    conda 'bioconda::vamb'

    publishDir path: "${params.out_dir}/vamb", mode: 'copy'

    input:
    tuple val(sample_id), path(polished_fasta), path(depth_file)

    output:
    tuple val(sample_id), path("vamb_bins/vamb_bin_*.fna"), emit: vamb_bins

    script:
    """
    # Ensure a clean output directory
    rm -rf vamb_bins && mkdir -p vamb_bins

    # Run Vamb with precomputed coverage and specified threads
    vamb bin default \
      --fasta ${polished_fasta} \
      --outdir vamb_bins \
      --minfasta 200000 \
      -p ${task.cpus}
    """
}


// -------------------------
// Process: runRosella (binning #4)
//    - Takes (sample_id, polished_fasta, depth_file) → emits (sample_id, rosella_bins/)
// -------------------------
process runRosella {
    tag { sample_id }
    cpus params.threads
    conda 'bioconda::rosella'

    publishDir path: "${params.out_dir}/rosella_primary", mode: 'copy'

    input:
    tuple val(sample_id), path(polished_fasta), path(depth_file)

    output:
    tuple val(sample_id), path("rosella_bins/rosella_bin_*.fna"), emit: rosella_bins

    script:
    """
    # Create output directory for bins
    mkdir -p rosella_bins

    # Run Rosella with precomputed coverage values (no CoverM step needed)
    rosella recover \
      -r ${polished_fasta} \
      --coverage-file ${depth_file} \
      -o rosella_bins \
      -t ${task.cpus}
    """
}


// -------------------------
// Process: runDASTool (ensemble binning)
//    - Takes (sample_id, semibinDir, metabatDir, vambDir, rosellaDir)
//      → emits (sample_id, dastool_bins/)
// -------------------------
process runDASTool {
    tag "${sample_id}"
    cpus params.threads
    conda 'bioconda::das_tool'

    publishDir path: "${params.out_dir}/dastool", mode: 'copy'

    input:
    tuple val(sample_id),
          path(semibinDir),
          path(metabatDir),
          path(vambDir),
          path(rosellaDir)

    output:
    tuple val(sample_id), path("dastool_bins"), emit: dastool_bins

    script:
    """
    mkdir -p dastool_bins
    
    # Create contig2bin files for each binner
    find ${semibinDir} -name "*.fa" -exec basename {} .fa \\; | awk -v OFS='\\t' '{print \$1, "semibin2"}' > semibin2.contigs2bin.tsv
    find ${metabatDir} -name "*.fa" -exec basename {} .fa \\; | awk -v OFS='\\t' '{print \$1, "metabat2"}' > metabat2.contigs2bin.tsv
    find ${vambDir} -name "*.fna" -exec basename {} .fna \\; | awk -v OFS='\\t' '{print \$1, "vamb"}' > vamb.contigs2bin.tsv
    find ${rosellaDir} -name "*.fna" -exec basename {} .fna \\; | awk -v OFS='\\t' '{print \$1, "rosella"}' > rosella.contigs2bin.tsv
    
    DAS_Tool \
      -i semibin2.contigs2bin.tsv,metabat2.contigs2bin.tsv,vamb.contigs2bin.tsv,rosella.contigs2bin.tsv \
      -l semibin2,metabat2,vamb,rosella \
      -c \$(find ${semibinDir} ${metabatDir} ${vambDir} ${rosellaDir} -name "*.fa" -o -name "*.fna" | head -1) \
      -o dastool_bins/DAS_Tool \
      --threads ${task.cpus}
    """
}


// -------------------------
// Process: runCheckM2 (QC on bins - final QC)
//    - Takes (sample_id, bins_to_qc) → emits (sample_id, checkm2_report.tsv)
// -------------------------
process runCheckM2 {
    tag "${sample_id}"
    cpus params.threads
    conda 'bioconda::checkm2'

    publishDir path: "${params.out_dir}/checkm2", mode: 'copy'

    input:
    tuple val(sample_id), path(bins_to_qc)

    output:
    tuple val(sample_id), path("checkm2_report.tsv"), emit: checkm2_tuple

    script:
    """
    checkm2 predict \
      --input ${bins_to_qc} \
      --output-directory checkm2_output \
      --threads ${task.cpus}
    
    cp checkm2_output/quality_report.tsv checkm2_report.tsv
    """
}

// -------------------------
// Process: runCheckM2_iter1 (QC on bins - iteration 1)
// -------------------------
process runCheckM2_iter1 {
    tag "${sample_id}"
    cpus params.threads
    conda 'bioconda::checkm2'

    publishDir path: "${params.out_dir}/checkm2_iter1", mode: 'copy'

    input:
    tuple val(sample_id), path(bins_to_qc)

    output:
    tuple val(sample_id), path("checkm2_report.tsv")

    script:
    """
    checkm2 predict \
      --input ${bins_to_qc} \
      --output-directory checkm2_output \
      --threads ${task.cpus}
    
    cp checkm2_output/quality_report.tsv checkm2_report.tsv
    """
}

// -------------------------
// Process: runCheckM2_iter2 (QC on bins - iteration 2)
// -------------------------
process runCheckM2_iter2 {
    tag "${sample_id}"
    cpus params.threads
    conda 'bioconda::checkm2'

    publishDir path: "${params.out_dir}/checkm2_iter2", mode: 'copy'

    input:
    tuple val(sample_id), path(bins_to_qc)

    output:
    tuple val(sample_id), path("checkm2_report.tsv")

    script:
    """
    checkm2 predict \
      --input ${bins_to_qc} \
      --output-directory checkm2_output \
      --threads ${task.cpus}
    
    cp checkm2_output/quality_report.tsv checkm2_report.tsv
    """
}

// -------------------------
// Process: runCheckM2_iter3 (QC on bins - iteration 3)
// -------------------------
process runCheckM2_iter3 {
    tag "${sample_id}"
    cpus params.threads
    conda 'bioconda::checkm2'

    publishDir path: "${params.out_dir}/checkm2_iter3", mode: 'copy'

    input:
    tuple val(sample_id), path(bins_to_qc)

    output:
    tuple val(sample_id), path("checkm2_report.tsv")

    script:
    """
    checkm2 predict \
      --input ${bins_to_qc} \
      --output-directory checkm2_output \
      --threads ${task.cpus}
    
    cp checkm2_output/quality_report.tsv checkm2_report.tsv
    """
}

// -------------------------
// Process: runCheckM2_iter4 (QC on bins - iteration 4)
// -------------------------
process runCheckM2_iter4 {
    tag "${sample_id}"
    cpus params.threads
    conda 'bioconda::checkm2'

    publishDir path: "${params.out_dir}/checkm2_iter4", mode: 'copy'

    input:
    tuple val(sample_id), path(bins_to_qc)

    output:
    tuple val(sample_id), path("checkm2_report.tsv")

    script:
    """
    checkm2 predict \
      --input ${bins_to_qc} \
      --output-directory checkm2_output \
      --threads ${task.cpus}
    
    cp checkm2_output/quality_report.tsv checkm2_report.tsv
    """
}

// -------------------------
// Process: runCheckM2_iter5 (QC on bins - iteration 5)
// -------------------------
process runCheckM2_iter5 {
    tag "${sample_id}"
    cpus params.threads
    conda 'bioconda::checkm2'

    publishDir path: "${params.out_dir}/checkm2_iter5", mode: 'copy'

    input:
    tuple val(sample_id), path(bins_to_qc)

    output:
    tuple val(sample_id), path("checkm2_report.tsv")

    script:
    """
    checkm2 predict \
      --input ${bins_to_qc} \
      --output-directory checkm2_output \
      --threads ${task.cpus}
    
    cp checkm2_output/quality_report.tsv checkm2_report.tsv
    """
}


// -------------------------
// Process: runRosellaRefine (refine bins - final refine)
//    - Takes (sample_id, checkm2_report.tsv) → emits (sample_id, refined_bins/)
// -------------------------
process runRosellaRefine {
    tag "${sample_id}"
    cpus params.threads
    conda 'bioconda::rosella'

    publishDir path: "${params.out_dir}/rosella_refine", mode: 'copy'

    input:
    tuple val(sample_id), path(checkm2_report)

    output:
    tuple val(sample_id), path("refined_bins"), emit: refined_bins

    script:
    """
    mkdir -p refined_bins
    rosella refine \
      --checkm-file ${checkm2_report} \
      --output refined_bins \
      --threads ${task.cpus}
    """
}

// -------------------------
// Process: runRosellaRefine_iter1 (refine bins - iteration 1)
// -------------------------
process runRosellaRefine_iter1 {
    tag "${sample_id}"
    cpus params.threads
    conda 'bioconda::rosella'

    publishDir path: "${params.out_dir}/rosella_refine_iter1", mode: 'copy'

    input:
    tuple val(sample_id), path(checkm2_report)

    output:
    tuple val(sample_id), path("refined_bins")

    script:
    """
    mkdir -p refined_bins
    rosella refine \
      --checkm-file ${checkm2_report} \
      --output refined_bins \
      --threads ${task.cpus}
    """
}

// -------------------------
// Process: runRosellaRefine_iter2 (refine bins - iteration 2)
// -------------------------
process runRosellaRefine_iter2 {
    tag "${sample_id}"
    cpus params.threads
    conda 'bioconda::rosella'

    publishDir path: "${params.out_dir}/rosella_refine_iter2", mode: 'copy'

    input:
    tuple val(sample_id), path(checkm2_report)

    output:
    tuple val(sample_id), path("refined_bins")

    script:
    """
    mkdir -p refined_bins
    rosella refine \
      --checkm-file ${checkm2_report} \
      --output refined_bins \
      --threads ${task.cpus}
    """
}

// -------------------------
// Process: runRosellaRefine_iter3 (refine bins - iteration 3)
// -------------------------
process runRosellaRefine_iter3 {
    tag "${sample_id}"
    cpus params.threads
    conda 'bioconda::rosella'

    publishDir path: "${params.out_dir}/rosella_refine_iter3", mode: 'copy'

    input:
    tuple val(sample_id), path(checkm2_report)

    output:
    tuple val(sample_id), path("refined_bins")

    script:
    """
    mkdir -p refined_bins
    rosella refine \
      --checkm-file ${checkm2_report} \
      --output refined_bins \
      --threads ${task.cpus}
    """
}

// -------------------------
// Process: runRosellaRefine_iter4 (refine bins - iteration 4)
// -------------------------
process runRosellaRefine_iter4 {
    tag "${sample_id}"
    cpus params.threads
    conda 'bioconda::rosella'

    publishDir path: "${params.out_dir}/rosella_refine_iter4", mode: 'copy'

    input:
    tuple val(sample_id), path(checkm2_report)

    output:
    tuple val(sample_id), path("refined_bins")

    script:
    """
    mkdir -p refined_bins
    rosella refine \
      --checkm-file ${checkm2_report} \
      --output refined_bins \
      --threads ${task.cpus}
    """
}

// -------------------------
// Process: runRosellaRefine_iter5 (refine bins - iteration 5)
// -------------------------
process runRosellaRefine_iter5 {
    tag "${sample_id}"
    cpus params.threads
    conda 'bioconda::rosella'

    publishDir path: "${params.out_dir}/rosella_refine_iter5", mode: 'copy'

    input:
    tuple val(sample_id), path(checkm2_report)

    output:
    tuple val(sample_id), path("refined_bins")

    script:
    """
    mkdir -p refined_bins
    rosella refine \
      --checkm-file ${checkm2_report} \
      --output refined_bins \
      --threads ${task.cpus}
    """
}


// -------------------------
// Process: runCoverMBins (coverage of refined bins vs. reads)
//    - Takes (sample_id, refined_bins, reads) → emits (sample_id, bins_depth.tsv)
// -------------------------
process runCoverMBins {
    tag "${sample_id}"
    cpus params.threads
    conda 'bioconda::coverm'

    publishDir path: "${params.out_dir}/coverm_bins", mode: 'copy'

    input:
    tuple val(sample_id), path(refined_bins), path(reads)

    output:
    tuple val(sample_id), path("bins_depth.tsv"), emit: final_depth_tuple

    script:
    def reads_list = (reads instanceof List) ? reads.join(' ') : reads.toString()
    """
    coverm genome \
      --coupled ${reads_list} \
      --genome-directory ${refined_bins} \
      --methods mean \
      --threads ${task.cpus} \
      --output-file bins_depth.tsv
    """
}


// -------------------------
// Process: runGtdbtk (classify bins)
//    - Takes (sample_id, refined_bins) → emits (sample_id, gtdbtk_report.tsv)
// -------------------------
process runGtdbtk {
    tag "${sample_id}"
    cpus params.threads
    conda 'bioconda::gtdbtk'

    publishDir path: "${params.out_dir}/gtdbtk", mode: 'copy'

    input:
    tuple val(sample_id), path(refined_bins)

    output:
    tuple val(sample_id), path("gtdbtk_report.tsv"), emit: gtdbtk_tuple

    script:
    """
    export GTDBTK_DATA_PATH=${params.gtdb_db}
    
    gtdbtk classify_wf \
      --genome_dir ${refined_bins} \
      --out_dir gtdbtk_output \
      --cpus ${task.cpus} \
      --extension fa

    # Combine bacterial and archaeal results if they exist
    if [ -f gtdbtk_output/gtdbtk.bac120.summary.tsv ]; then
        cp gtdbtk_output/gtdbtk.bac120.summary.tsv gtdbtk_report.tsv
    elif [ -f gtdbtk_output/gtdbtk.ar53.summary.tsv ]; then
        cp gtdbtk_output/gtdbtk.ar53.summary.tsv gtdbtk_report.tsv
    else
        touch gtdbtk_report.tsv
    fi
    """
}


// -------------------------
// Process: runMultiQC (aggregate results)
//    - Takes a single 6‐element tuple:
//      (sample_id, singlem_report, depth_file, checkm2_report, final_depth, gtdbtk_report)
//    - Emits (sample_id, multiqc_report.html)
// -------------------------
process runMultiQC {
    tag "${sample_id}"
    cpus params.threads
    conda 'bioconda::multiqc'

    publishDir path: "${params.out_dir}/final_report", mode: 'copy'

    input:
    tuple val(sample_id),
          path(singlem_report),
          path(depth_file),
          path(checkm2_report),
          path(final_depth),
          path(gtdbtk_report)

    output:
    tuple val(sample_id), path("multiqc_report.html"), emit: multiqc_tuple

    script:
    """
    # Create a directory with all input files for MultiQC
    mkdir -p multiqc_input
    cp ${singlem_report} multiqc_input/
    cp ${depth_file} multiqc_input/
    cp ${checkm2_report} multiqc_input/
    cp ${final_depth} multiqc_input/
    cp ${gtdbtk_report} multiqc_input/
    
    multiqc multiqc_input -o . -n multiqc_report.html
    """
}

// -------------------------
// Sub‐workflow: refine_qc (5 iterations of CheckM2 + Rosella refine)
// -------------------------
workflow refine_qc {
    take:
    bins_ch   // expects (sample_id, binsDir)

    main:
    // First iteration
    checkm2_1 = runCheckM2_iter1(bins_ch)
    refined_1 = runRosellaRefine_iter1(checkm2_1)
    
    // Second iteration
    checkm2_2 = runCheckM2_iter2(refined_1)
    refined_2 = runRosellaRefine_iter2(checkm2_2)
    
    // Third iteration
    checkm2_3 = runCheckM2_iter3(refined_2)
    refined_3 = runRosellaRefine_iter3(checkm2_3)
    
    // Fourth iteration
    checkm2_4 = runCheckM2_iter4(refined_3)
    refined_4 = runRosellaRefine_iter4(checkm2_4)
    
    // Fifth iteration
    checkm2_5 = runCheckM2_iter5(refined_4)
    final_refined = runRosellaRefine_iter5(checkm2_5)

    emit:
    refined_bins = final_refined  // emits (sample_id, refined_binsDir)
}

workflow {
    // Normalize out_dir to avoid redefinition warning
    def normalized_out_dir = params.out_dir.replaceAll('/+$','')
    
    // 1) Create channel from folder with *.fastq.gz files
    // Option 1: If you have paired-end reads with _1/_2 or _R1/_R2 naming
    if (params.containsKey('paired') && params.paired) {
        reads_ch = Channel
            .fromFilePairs("${params.input_dir}/*_{1,2}.fastq.gz")
            .ifEmpty { error "No paired-end files found in ${params.input_dir} with pattern *_{1,2}.fastq.gz" }
            .map { sample_id, files ->
                tuple(sample_id, files)
            }
    } else {
        // Option 2: Single sample with all fastq.gz files (for single-end or mixed)
        reads_ch = Channel
            .fromPath("${params.input_dir}/*.fastq.gz")
            .ifEmpty { error "No fastq.gz files found in ${params.input_dir}" }
            .collect()
            .map { files ->
                def sample_id = file(params.input_dir).name
                tuple(sample_id, files)
            }
    }

    // Debug: Print what files are found
    reads_ch.view { sample_id, files ->
        "Found sample: ${sample_id} with ${files.size()} files: ${files}"
    }

    // 2) Run MyLOASM → emits polished_tuple = (sample_id, assembly_primary.fa)
    polished_tuple = runMyloasm(reads_ch)

    // 3) Run SingleM → emits singlem_tuple = (sample_id, singlem_report.tsv)
    singlem_tuple = runSingleM(polished_tuple)

    // 4) Run CoverM by joining polished_tuple and reads_ch → emits (sample_id, depth.tsv)
    coverm_input_ch = polished_tuple
        .join(reads_ch)
        .map { sample_id, faFile, readsFiles ->
            tuple(sample_id, faFile, readsFiles)
        }
    depth_tuple = runCoverM(coverm_input_ch)

    // 5) Prepare input for all binning tools: (sample_id, polished_fasta, depth_file)
    semibin_input_ch = polished_tuple
        .join(depth_tuple)
        .map { sample_id, faFile, depthFile ->
            tuple(sample_id, faFile, depthFile)
        }

    // 6) Run each binning tool in parallel:
    semibin_bins   = runSemibin2(semibin_input_ch)
    metabat2_bins  = runMetabat2(semibin_input_ch)
    vamb_bins      = runVamb(semibin_input_ch)
    rosella_bins   = runRosella(semibin_input_ch)

    // 7) Join the four binners by sample_id:
    //    (sample_id, semibinDir, metabatDir, vambDir, rosellaDir)
    das_input_ch = semibin_bins
        .join(metabat2_bins)
        .join(vamb_bins)
        .join(rosella_bins)
        .map { sample_id, semibinDir, metabatDir, vambDir, rosellaDir ->
            tuple(sample_id, semibinDir, metabatDir, vambDir, rosellaDir)
        }
    dastool_bins = runDASTool(das_input_ch)

    // 8) Sub‐workflow: iterative QC & refine (5 rounds)
    refined_bins = refine_qc(dastool_bins)

    // 9) One final QC & coverage on refined bins:
    checkm2_tuple = runCheckM2(refined_bins)
    coverm_bins_input_ch = refined_bins
        .join(reads_ch)
        .map { sample_id, binsDir, readsFiles ->
            tuple(sample_id, binsDir, readsFiles)
        }
    final_depth_tuple = runCoverMBins(coverm_bins_input_ch)

    // 10) Classify refined bins with GTDB-Tk
    gtdbtk_tuple = runGtdbtk(refined_bins)

    // 11) Aggregate all results with MultiQC:
    //     Channel must emit (sample_id, singlem_report, depth, checkm2_report, final_depth, gtdbtk_report)
    multiqc_input_ch = singlem_tuple
        .join(depth_tuple)
        .join(checkm2_tuple)
        .join(final_depth_tuple)
        .join(gtdbtk_tuple)
        .map { sample_id, singleReport, depthReport, checkm2Report, finalDepthReport, gtdbtkReport ->
            tuple(sample_id, singleReport, depthReport, checkm2Report, finalDepthReport, gtdbtkReport)
        }
    multiqc_tuple = runMultiQC(multiqc_input_ch)

    // 12) Print out final MultiQC report for each sample
    multiqc_tuple.view { sample_id, report ->
        "Sample ${sample_id} → Final report: ${report}"
    }
}
