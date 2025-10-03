# Somatem
LLM accessible long-read metagenomics pipeline with best practices
## Outline
[Plannning tools: for development](https://github.com/ppreshant/somatem-docs/tree/main?tab=readme-ov-file#plannning-tools) | [Overarching goals](https://github.com/treangenlab/SOMAteM/#overarching-goals)

Somatem is a collection of nextflow scripts for long-read sequencing metagenomics from 4th generation sequencing technologies (Oxford Nanopore Technologies and PacBio). The scripts are designed to be modular and easy to use, with a focus on best practices for long-read sequencing data analysis.

The pipeline includes key subworkflows for 
- data pre-processing: Visualize quality with [NanoPlot](https://github.com/wdecoster/NanoPlot), then remove host contamination using [hostile](https://github.com/bede/hostile), followed by sequence quality and length filtering with [chopper](https://github.com/wdecoster/chopper), then one last __NanoPlot__ analysis to show the statistics of your final reads
- taxonomic classification: 
- genome assembly / metagenome-assembled genome (MAG) analysis: 
- pathogen detection: using seqscreen to screen the reads for pathogens: 



## Initial setup

1. Clone this repo from GitHub
```
git clone https://github.com/treangenlab/Somatem
```

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
cp docs/somatem_docs/metadata_template.yaml assets/custom_metadata.yaml
```

- Run the pipeline from the base directory (`Somatem/`) using the following command in the terminal:
```bash
nextflow run . -param-file assets/custom_metadata.yaml
```

Note that the `assembly_mags` section could take a few hours to run and it depends on the complexity of your data and your computational resources. It took 6 hours to run the 2 example files `assets/mag_big_samplesheet.csv` on a cluster with 128 GB memory, 128 cpus and 6 TB of free space.

- This command automatically downloads required databases ranging from 2-60 GB size: 
For the `assembly_mags` workflow there are a few large sized databases including the [checkm2](https://github.com/chklovski/CheckM2) and [gtdbtk](https://gtdb.ecogenomic.org/) databases, as well as the [singlem metapackage](https://zenodo.org/records/15232972).

**Cleanup the README below this**

GTDB-Tk v2.4.1 requires ~140G of external data which needs to be downloaded and extracted. This can be done automatically, or manually.
```
# Automatic:
# Run the command "download-db.sh" to automatically download and extract to:
/path/to/miniforge3/envs/singlem/share/gtdbtk-2.4.1/db/

# Manual:
# Manually download the latest reference data:
wget https://data.ace.uq.edu.au/public/gtdb/data/releases/release226/226.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r226_data.tar.gz

# Extract the archive to a target directory:
tar -xvzf gtdbtk_r226_data.tar.gz -C "/path/to/target/db" --strip 1 › /dev/null
rm gtdbtk_r226_data.tar.gz

# Set the GTDBTK DATA PATH environment variable by running:
conda env config vars set GTDBTK DATA PATH="/path/to/target/db"
```


---
# Plannning and other documentation
_We will make a wiki to document the planning tools to be included in the pipeline ([flowchart](../docs/SOMAteM-sketch-v1.2.jpg)). This can be moved to the [wiki](https://docs.github.com/en/communities/documenting-your-project-with-wikis/adding-or-editing-wiki-pages#cloning-wikis-to-your-computer) when it is created eventually._

See `docs/` for more documentation stuff including: [Plannning tools](https://github.com/treangenlab/SOMAteM/blob/main/docs/bioinformatic_tools_planner.md#plannning-tools)

---

# Overarching goals

**Outline**: Goal is to make a novel bioinformatic mega-pipeline that incorporates best practices the key contribution of the paper. We will pivot omi's current implementation (RAG-LLM) to work as a UI to run this pipeline. Omi will use text prompts to choose between various paths of tool-calls within the flowchart of this nextflow pipeline. The pipeline will have some decision points to choose between overarching themes (such as `genome assembly` vs using `reads` directly) and between different tools with proficiencies in speed, accuracy, and false-positives. 

**Points:**
- **Novelty/unfilled niche**: The pipeline will be geared towards the newer long-read technologies (`nanopore`, `pacbio`) where there aren't many established/optimal tools yet
- **Better bioinformatic tools**: Todd's lab has a lot of tools that are fast and efficient for this. 
        - This project will be a good wrapper for promoting all the good tools produced by the lab ; and benefits of LLM
- Clear goal & **timeline**: A concrete pipeline allows us to move faster towards a paper (*outline around Sep 1st*)
 - **Competition**: Having our own custom made pipeline gives  an edge over the [highly resourced](https://seqera.io/blog/seqera-raises-26m-series-b/) seqera company's (owns nextflow) [AI-chat](https://seqera.io/ask-ai/chat) doing [similar work](https://www.healthcareittoday.com/2024/08/29/seqera-acquires-tinybio-to-advance-science-for-everyone-now-through-genai/) (news [release](https://seqera.io/blog/seqera-ai-launch/))
