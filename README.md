# somatem
**A modular and open source metagenomic analysis toolkit designed for long reads**

somatem is a modular Nextflow based pipeline designed for long-read microbiome analysis, including both 16S and metagenomic support. somatem supports both Oxford Nanopore Technologies and PacBio. Built with ease of use and analytical rigor in mind, somatem enforces best practices for long-read sequencing data analysis.

The pipeline is divided into key subworkflows, allowing users to run the exact analyses they need:
* **Pre-processing:** Quality control and read filtering.
* **Taxonomic Profiling:** Taxonomic classification and relative abundance estimation.
* **Assembly & MAG Analysis:** *De novo* metagenomic assembly, binning, quality assessment, and functional annotation.
* **Genome Dynamics:** Structural variant and horizontal gene transfer detection for temporal samples.

---

## Initial Setup

Follow these steps to configure your environment and download the somatem pipeline. Note: This pipeline is designed for Linux/macOS environments and is not compatible with Windows.

**1. Clone the Repository**
Clone the somatem repository along with its required submodules:
```bash
git clone --recurse-submodules https://github.com/treangenlab/somatem
cd somatem
```
> **Note:** If you accidentally cloned the repository without the `--recurse-submodules` flag, you can fetch them by running `git submodule update --init --recursive` inside the repo directory. See the [Git Submodules documentation](https://git-scm.com/book/en/v2/Git-Tools-Submodules) for troubleshooting.

**2. Install Micromamba**
We utilize `micromamba` (a faster, drop-in replacement for `conda`) to manage environments. Install the latest version:
```bash
"${SHELL}" <(curl -L https://micro.mamba.pm/install.sh)
```

**3. Create and Activate the Nextflow Environment**
Set up a dedicated base environment for Nextflow: (this installs nextflow and nf-core; _needs to be run once_)
```bash
micromamba create -f nf_base_env.yml
```
Activate the environment each time you use the pipeline:
```bash
micromamba activate nf_base_env
```

**4. Download Example Data (Optional but Recommended)**
To verify your installation, you can download our provided test datasets. From the root `somatem/` directory (with your environment activated), run:
```bash
nextflow run subworkflows/local/get_example_data.nf
```
*This will populate the `examples/data` directory with sample files for pipeline testing.*

---

## Database Configuration

Several tools in this pipeline rely on large reference databases. Proper configuration is essential to manage storage effectively.

* **Storage Requirements:** Some databases (e.g., Bakta, CheckM2, SingleM) require up to 100 GB of free space. Ensure your target drive has adequate capacity.
* **Directory Setup:** Decide whether you want a local database directory within the `somatem` folder, or a shared, centralized directory (highly recommended for HPC cluster environments).
* **Configuration:** Update the `nextflow.config` file to point the pipeline to your chosen directory. Locate and modify the following variable:
    ```groovy
    db_base_dir = "/home/dbs" // Change this to "./assets/databases" for local storage
    ```

---

## Usage

Due to the modular design of somatem, you must configure the pipeline to run the specific subworkflows relevant to your research questions.

**1. Prepare Your Metadata**
Copy the provided metadata template to create your custom configuration file. The template contains detailed comments to guide you through the available parameters.
```bash
cp docs/somatem-docs/metadata_template.yaml assets/custom_metadata.yaml
```

**2. Execute the Pipeline**
Run somatem from the base directory, passing in your customized metadata file:
```bash
nextflow run . -param-file assets/custom_metadata.yaml
```

**Performance & Resource Notes:**
* **Automated Downloads:** The pipeline automatically downloads most required databases (<3 GB). However, the [Bakta database](https://zenodo.org/records/14916843) used in the `assembly_mags` subworkflow is approximately 60 GB and may require additional time.
* **Compute Time:** The `assembly_mags` step is computationally intensive. As a benchmark, processing the two example files (`assets/mag_big_samplesheet.csv`) takes roughly 6 hours on an HPC cluster equipped with 128 CPUs, 128 GB of memory, and 2 TB of free storage.

---

## Pipeline Tools

somatem integrates state-of-the-art bioinformatics tools, neatly organized into the following subworkflows:

### Pre-processing
Prepares raw data for downstream analysis through rigorous quality control and filtering.
* **[NanoPlot](https://github.com/wdecoster/NanoPlot):** QC plotting suite for initial and final assessment of long-read sequencing data.
* **[Hostile](https://github.com/bede/hostile):** Depletes host contamination by filtering reads that align to a host reference genome.
* **[Chopper](https://github.com/wdecoster/chopper):** Filters nanopore reads by quality and length, removing sub-par data.

### Taxonomic Profiling
Delivers rapid and accurate taxonomic classification for metagenomic datasets.
* **[Emu](https://github.com/treangenlab/emu):** Taxonomic classification and abundance estimation optimized for long-read 16S rRNA.
* **[Lemur](https://github.com/treangenlab/lemur):** Rapid, multi-marker gene taxonomic profiling for long-read metagenomes.
* **[MAGnet](https://github.com/treangenlab/magnet):** Refines taxonomic profiles via reference genome mapping to correct false positives.
* **[SingleM](https://github.com/wwood/singlem):** Profiles microbial communities using universal marker genes. Includes the `pipe` module for reads/assemblies and the `appraise` module to evaluate binning completeness.

### Assembly & MAG Analysis
Handles *de novo* assembly, genome binning, and functional annotation.
* **[Flye](https://github.com/fenderglass/Flye):** Repeat-graph-based *de novo* assembler optimized for PacBio and Nanopore reads.
* **[Minimap2](https://github.com/lh3/minimap2) & [SAMtools](https://github.com/samtools/samtools):** Pairwise alignment processing, read mapping, and coverage calculation.
* **[SemiBin2](https://github.com/BigDataBiology/SemiBin):** Metagenomic binning leveraging semi-supervised deep learning.
* **[CheckM2](https://github.com/chklovski/CheckM2):** Machine-learning-driven prediction of genome bin quality and completeness.
* **[Bakta](https://github.com/oschwengers/bakta):** Comprehensive and rapid annotation of bacterial genomes and plasmids.

### Genome Dynamics
Investigates structural variations over time.
* **[Rhea](https://github.com/treangenlab/rhea):** Detects structural variants and horizontal gene transfer events in temporally evolving microbial samples.
* **[Bandage](https://github.com/rrwick/Bandage):** Interactive visualization tool for assembly graphs, highly useful for reviewing Rhea outputs.

### Functional Annotation
Screens for targets of clinical and functional interest.
* **[SeqScreen](https://gitlab.com/treangenlab/seqscreen):** Functional screening of pathogenic sequences and antimicrobial resistance (AMR) genes.

### Reporting & Visualization
Aggregates and visualizes complex datasets.
* **[Taxburst](https://github.com/taxburst/taxburst):** Interactive, web-based visualization of taxonomic profiles.
* **[MultiQC](https://github.com/multiqc/multiqc):** Aggregates logs and results across multiple tools into a single, user-friendly HTML report.

---

## Additional Documentation

For deeper dives into pipeline architecture and tool notes, please see the `docs/` directory:
* [Installation Guide](docs/installation.md)
* [Planning Tools & Development Roadmap: ARCHIVED](https://github.com/treangenlab/somatem-docs/blob/main/planning/bioinformatic_tools_planner.md)
* [Tool Status and Quick Notes](docs/tool_status-quicknotes.md)

## Citation

If somatem facilitates your research, please cite the underlying tools that made your analysis possible. A comprehensive list of citation links is available in [docs/somatem-docs/tool_links.csv](https://github.com/treangenlab/somatem-docs/blob/main/tool_links.csv).

## Contributing & License

Contributions from the community are welcome! Please review our development documentation for guidelines on how to submit pull requests. 

This project is licensed under the **GNU General Public License v3.0 (GPLv3)**. See the `LICENSE` file for full details.
