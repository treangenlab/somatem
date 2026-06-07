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

### 1. Install conda/mamba/micromamba

We utilize `micromamba` (a faster, drop-in replacement for `conda`) but any of the listed package managers will work for to install somatem. Install micromamba using the command below in Linux. Source: [docs](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html)

```sh
"${SHELL}" <(curl -L micro.mamba.pm/install.sh)
```

### 2. Create and Activate the somatem Environment
Set up a dedicated base environment for somatem:

```bash
micromamba create -n somatem -c bioconda somatem # Again use your package manager of interest
```

### 3. Set environment variables (optional. only if changing defaults)
You may want to override the default locations for somatem's database storage, conda environment cache, or other data.

You can set the following environment variables:
- `SOMATEM_DB_DIR`: Directory for downloaded databases (default: `~/somatem_databases`). _Change this if you want to store databases in a shared location with other users/other projects etc. to minimize storage and if running on HPC clusters_
- `NXF_CONDA_CACHEDIR`: Directory for conda environment cache (default: `~/.nextflow/cache`). _Change this if you want to reallocate storage into scratch or something if on a HPC cluster. Note that When using a computing cluster it must be a shared folder accessible from all compute nodes._
- `SOMATEM_UNIFIED_DB_DIR`: Directory for unified database files for ensemble species detection (default: same as `SOMATEM_DB_DIR`). _note: this is a temporary location. These DBs will eventually be integrated into the db dir and this variable will be removed_

Environment variables can be set by exporting using `export SOMATEM_DB_DIR=/path/to/dbs` in the terminal.
You can edit the paths in `assets/scripts/somatem_env.sh` and add it your shell's profile file (e.g., `.bashrc`, `.zshrc`) so it's loaded automatically for future logins:

```bash
echo "source /path/to/somatem_env.sh" >> ~/.bashrc
```

### 4. Test out somatem!

To process long-read 16S sequencing with somatem one would simply
Activate the environment using this before running somatem: (run on each new terminal session)

```bash
# activate environment
micromamba activate somatem

# run the somatem 16S subworkflow
somatem 16S -i /path/to/16S_samplesheet.csv -o /path/to/desired_output
```
For help on making your input samplesheet, please see the example [here](https://github.com/treangenlab/somatem/blob/main/assets/16S_sheet.csv)

Note: if you are actively developing the pipeline, check out docs/dev-notes.md(docs/dev-notes.md) for extra setup instructions.
---

## Usage

Information on how to run the various subworkflows in somatem can be found in our [wiki pages](https://github.com/treangenlab/somatem/wiki)!

## Database Configuration

Several tools in this pipeline rely on large reference databases. Proper configuration is essential to manage storage effectively. The first time you run a pipeline requiring a database these will be installed for you and saved at that path for future runs.

* **Storage Requirements:** Some databases (e.g., Bakta, CheckM2, SingleM) require up to 100 GB of free space. Ensure your target drive has adequate capacity.
* **Directory Setup:** By default, Somatem stores generated Nextflow conda environments and downloaded databases under the active conda environment at `$CONDA_PREFIX/var/somatem`. If no conda environment is active, it falls back to `~/.somatem`.
* **Configuration:** To use another location, set one of these environment variables before running Somatem:
    ```bash
    export SOMATEM_HOME=/path/to/somatem-data      # sets both databases and conda cache
    export SOMATEM_DB_DIR=/path/to/databases       # overrides databases only
    export SOMATEM_CONDA_CACHE=/path/to/nxf-envs   # overrides Nextflow conda environments only
    export SOMATEM_UNIFIED_DB_DIR=/path/to/unified # overrides ensemble profiling databases
    ```

---

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
* [Planning Tools & Development Roadmap: ARCHIVED](https://github.com/treangenlab/somatem-docs/blob/main/planning/bioinformatic_tools_planner.md)
* [Tool Status and Quick Notes](docs/tool_status-quicknotes.md)

## Citation

If somatem facilitates your research, please cite the underlying tools that made your analysis possible. A comprehensive list of citation links is available in [docs/somatem-docs/tool_links.csv](https://github.com/treangenlab/somatem-docs/blob/main/tool_links.csv).

## Contributing & License

Contributions from the community are welcome! Please review our development documentation for guidelines on how to submit pull requests. 

This project is licensed under the **GNU General Public License v3.0 (GPLv3)**. See the `LICENSE` file for full details.
