# Somatem
LLM accessible long-read metagenomics pipeline with best practices
## Outline
[Plannning tools](https://github.com/ppreshant/somatem-docs/tree/main?tab=readme-ov-file#plannning-tools) | [Overarching goals](https://github.com/treangenlab/SOMAteM/#overarching-goals)

## Installation

To run the workflows you will need to clone the github and make conda environments for each of the workflows.

For the `somatem_prep.nf` workflow, the install is as follows:
```
# have a working conda install (module load conda or download sh file)

# clone this branch of the repo and move to dir 
git clone https://github.com/treangenlab/SOMAteM -b agm && cd SOMAteM

# create the conda env with all the tools
conda env create -n somatem_prep -f envs/somatem_prep.yml

# activate the conda env
conda activate somatem_prep

# and you should be good to run!
```

## Usage

#### long read preprocessing pipeline
We are currently developing nextflow scripts for commons long-read sequencing tasks. The first of which is called `somatem_prep.nf` which will take in fastq files (found in examples/data) and run [NanoPlot](https://github.com/wdecoster/NanoPlot) on the raw files, then remove host contamination using [hostile](https://github.com/bede/hostile), followed by sequence quality and length filtering with [chopper](https://github.com/wdecoster/chopper), then one last __NanoPlot__ analysis to show the statistics of your final reads, ready for use in taxonomic classification and assembly pipelines _under active development_.

`somatem_prep.nf`
```
# activate conda env
conda activate somatem_prep

# test on our example data
gzip -d examples/data/*.fastq.gz

# for simplicity one can use a premade config file
nextflow run /path/to/SOMAteM/workflows/somatem_prep.nf -c /path/to/SOMAteM/confs/somatem_prep.config
```
#### long read metagenome assembled genome pipeline
A custom nextflow pipeline of the current best practice tools for long read MAG assembly and binning. This pipeline is much more compute intensive. 

`somatem_mags.nf`
```
# assuming you already have nextflow and conda installed along with git repo cloned

# go to SOMAteM repo
cd SOMAteM 

# lil test run
nextflow run workflows/somatem_mags.nf --input_dir /path/to/SOMAteM/examples/data --output_dir /path/to/SOMAteM/examples/soma_mags_test --threads 8 --gtdb_db /path/to/gtdb214 -c /path/to/SOMAteM/confs/somatem_mags.config
```

For the `somatem_mags.nf` workflow there are also some databases you must download before use including the [checkm2](https://github.com/chklovski/CheckM2) and [gtdbtk](https://gtdb.ecogenomic.org/) databases, as well as the [singlem metapackage](https://zenodo.org/records/15232972).

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

### the somatem_mags.nf workflow still has some bugs but are being worked out in the very near future!!!

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
