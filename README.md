# Somatem
LLM accessible long-read metagenomics pipeline with best practices
## Outline
[Planning tools: for development](https://github.com/ppreshant/somatem-docs/tree/main?tab=readme-ov-file#planning-tools) | [Overarching goals](https://github.com/treangenlab/SOMAteM/#overarching-goals)

Somatem is a collection of nextflow scripts for long-read sequencing metagenomics from 4th generation sequencing technologies (Oxford Nanopore Technologies and PacBio). The scripts are designed to be modular and easy to use, with a focus on best practices for long-read sequencing data analysis.

The pipeline includes key [subworkflows](https://github.com/treangenlab/SOMAteM/tree/main/?tab=readme-ov-file#pipeline-tools) elaborated below for 
- **Pre-processing**: Quality control and read filtering
- **Taxonomic profiling**: Classification and abundance estimation
- **Assembly & MAG analysis**: De novo assembly, binning, quality assessment, and functional annotation
- **Genome dynamics**: Structural variant and horizontal gene transfer detection for timecourse samples


## Initial setup

1. Clone this repo from GitHub
```
git clone --recurse-submodules https://github.com/treangenlab/Somatem
```
    - Note: The repo contains submodules that need to be cloned as well. If you cloned without the `--recurse-submodules` flag, then run `git submodule update --init --recursive` in the repo directory. Read [submodules documentation](https://git-scm.com/book/en/v2/Git-Tools-Submodules) for more information or troubleshooting.

2. Install the latest version of micromamba from [here](https://github.com/mamba-org/mamba/releases) using `"${SHELL}" <(curl -L https://micro.mamba.pm/install.sh)` for linux/MacOS (This pipeline cannot run on Windows). `micromamba` is a faster version of `conda` that is used to create and manage conda environments. _We will use conda/micromamba terms interchangably_

3. Navigate to the cloned repo directory `Somatem/` and create a conda environment for nextflow. 
```
micromamba create -n nf_base_env nextflow
```

4. Activate the conda environment and install nextflow
```
micromamba activate nf_base_env
micromamba install nextflow
```
Now you are ready to run nextflow. If you want to test the pipelines with example datasets, proceed to step 5.

5. Download the example data using the `get_example_data` subworkflow by running this command in the root directory of the repo (after activating the `nf_base_env` conda environment):
```bash
micromamba activate nf_base_env
nextflow run subworkflows/local/get_example_data.nf
```
_This will download the example data to the `examples/data` directory for testing the pipelines._

## Database configuration
1. Decide if you want a dedicated directory for the databases within the `Somatem` directory or use a common shared directory for other members of your organization (highly recommended if you are using a cluster). 

2. Some of the databases (e.g. bakta, checkm2, singlem) run upto 100GB of space, so make sure you have enough space available

3. Update the `nextflow.config` file to point to the database directory by modifying this variable
```
db_base_dir                = "/home/dbs"
```
Change this to `./assets/databases` if not using a shared directory

## Usage
Since the pipeline has multiple subworkflows, you have to pick which one(s) are most relevant to your use case. 

- Input your configuration by copying the `metadata_template.yaml` file and filling it out. _The file has helpful comments to guide you through the process_
```sh
cp docs/somatem-docs/metadata_template.yaml assets/custom_metadata.yaml
```

- Run the pipeline from the base directory (`Somatem/`) using the following command in the terminal:
```bash
nextflow run . -param-file assets/custom_metadata.yaml
```

Note that the `assembly_mags` section could take a few hours to run and it depends on the complexity of your data and your computational resources. It took 6 hours to run the 2 example files `assets/mag_big_samplesheet.csv` on a cluster with 128 GB memory, 128 cpus and 6 TB of free space.

- This command automatically downloads required databases ranging from <3 GB size, except [bakta](https://zenodo.org/records/14916843) used in `assembly_mags` which is a large database (60GB of space).

## Pipeline Tools

Somatem integrates state-of-the-art bioinformatics tools organized into modular subworkflows:

### Pre-processing
Quality control and read filtering to prepare data for downstream analysis.

- **[NanoPlot](https://github.com/wdecoster/NanoPlot)** - QC plotting suite for long-read sequencing data and alignments. Used for initial and final quality assessment.
- **[hostile](https://github.com/bede/hostile)** - Removes host contamination from microbial metagenomes by filtering reads that align to a host genome.
- **[chopper](https://github.com/wdecoster/chopper)** - Filters nanopore sequencing reads by quality and length to remove low-quality and short reads.

### Taxonomic Profiling
Rapid and accurate taxonomic classification for both 16S amplicon and metagenomic data.

- **[Emu](https://github.com/treangenlab/emu)** - Taxonomic classification and abundance estimation for 16S rRNA reads optimized for long-read data.
- **[Lemur](https://github.com/treangenlab/lemur)** - Rapid and accurate taxonomic profiling on long-read metagenomic datasets using multi-marker genes.
- **[MAGnet](https://github.com/treangenlab/magnet)** - Refines taxonomic profiles for accuracy using reference genome mapping to correct false positives.
- **[SingleM pipe](https://github.com/wwood/singlem)** - Profiles microbial communities using universal marker genes for both reads and assembled genomes.
- **[SingleM appraise](https://github.com/wwood/singlem)** - Assesses and compares metagenomic samples to evaluate binning completeness.

### Assembly & MAG Analysis
De novo assembly, binning, quality assessment, and functional annotation.

- **[Flye](https://github.com/fenderglass/Flye)** - De novo assembler for single-molecule sequencing reads using repeat graphs. Optimized for PacBio and Oxford Nanopore.
- **[minimap2](https://github.com/lh3/minimap2)** - Pairwise alignment of assemblies to reference genomes and read mapping.
- **[samtools](https://github.com/samtools/samtools)** - Alignment processing and coverage calculation.
- **[SemiBin2](https://github.com/BigDataBiology/SemiBin)** - Metagenomic binning with semi-supervised deep learning for improved binning accuracy.
- **[CheckM2](https://github.com/chklovski/CheckM2)** - Accurate prediction of genome quality using machine learning for fast bin quality assessment.
- **[bakta](https://github.com/oschwengers/bakta)** - Rapid and accurate annotation of bacterial genomes and plasmids. Comprehensive annotation pipeline for high-quality bins.

### Genome Dynamics
Detection of structural variants and horizontal gene transfer events.

- **[rhea](https://github.com/treangenlab/rhea)** - Detects structural variants and horizontal gene transfer between temporally evolving microbial metagenomic samples.
- **[bandage](https://github.com/rrwick/Bandage)** - Interactive visualization of assembly graphs. Useful for visualizing results after rhea analysis.

### Functional Annotation
Screening for pathogenic sequences and antimicrobial resistance genes.

- **[SeqScreen](https://gitlab.com/treangenlab/seqscreen)** - Functional screening of pathogenic sequences in metagenomic data, including antibiotic resistance genes.

### Reporting & Visualization
Aggregation and interactive visualization of analysis results.

- **[taxburst](https://github.com/taxburst/taxburst)** - Interactive web-based visualization of taxonomic profiles.
- **[MultiQC](https://github.com/multiqc/multiqc)** - Aggregates and visualizes results from multiple tools and samples in a single report.


---

## Additional Documentation

For more detailed information, see the `docs/` directory:
- [Planning tools and development roadmap](https://github.com/treangenlab/SOMAteM/blob/main/docs/bioinformatic_tools_planner.md)
- [Installation guide](docs/installation.md)
- [Tool status and notes](docs/tool_status-quicknotes.md)

## Citation

If you use Somatem in your research, please cite the tools used in your analysis. Links to citations are available in the [tool_links.csv](docs/somatem-docs/tool_links.csv) file.

## Contributing

Contributions are welcome! Please see our development documentation for guidelines.

## License

This project is licensed under the GPL-3.0 License - see the LICENSE file for details.
