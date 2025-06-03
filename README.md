# SOMAteM
LLM accessible long-read metagenomics pipeline with best practices
## Outline
[Plannning tools](https://github.com/treangenlab/SOMAteM/#plannning-tools) | [Overarching goals](https://github.com/treangenlab/SOMAteM/#overarching-goals)


## Usage thus far

We are currently developing nextflow scripts for commons long-read sequencing tasks. The first of which is called `somatem_prep.nf` which will take in fastq files (found in examples/data) and run [NanoPlot](https://github.com/wdecoster/NanoPlot) on the raw files, then remove host contamination using [hostile](https://github.com/bede/hostile), followed by sequence quality and length filtering with [chopper](https://github.com/wdecoster/chopper), then one last __NanoPlot__ analysis to show the statistics of your final reads, ready for use in taxonomic classification and assembly pipelines _under active development_.

### Usage:

`somatem_prep.nf`
```
# to specify your own parameters within the CLI you can run something like this
nextflow run /path/to/SOMAteM/workflows/somatem_prep.nf --input_dir /path/to/SOMAteM/examples/data --output_dir /path/to/SOMAteM/examples/soma_prep_out --threads 12 --maxlength 30000 --minq 10 --minlen 250 --host_index 'human-t2t-hla'

# or for simplicity one can use a premade config file
nextflow run /path/to/SOMAteM/workflows/somatem_prep.nf -c /path/to/SOMAteM/confs/somatem_prep.config
```
*more to come soon!!!*

# Plannning tools
_We will make a wiki to document the planning tools to be included in the pipeline ([flowchart](../docs/SOMAteM-sketch-v1.2.jpg)). This can be moved to the [wiki](https://docs.github.com/en/communities/documenting-your-project-with-wikis/adding-or-editing-wiki-pages#cloning-wikis-to-your-computer) when it is created eventually._

## Major decision points
- 1 vs n samples
  - Rheaa is good for n samples so it will dictate the rest of the pipeline 
- Assembly vs reads
- Reference guided vs de novo assembly
- 3 broad branches of the pipeline
  - 16S
  - metagenomic classification
  - assembly based metagenomics


---

# Overarching goals

**Outline**: Goal is to make a novel bioinformatic mega-pipeline that incorporates best practices the key contribution of the paper. We will pivot omi's current implementation (RAG-LLM) to work as a UI to run this pipeline. Omi will use text prompts to choose between various paths of tool-calls within the flowchart of this nextflow pipeline. The pipeline will have some decision points to choose between overarching themes (such as `genome assembly` vs using `reads` directly) and between different tools with proficiencies in speed, accuracy, and false-positives. 

**Points:**
- **Novelty/unfilled niche**: The pipeline will be geared towards the newer long-read technologies (`nanopore`, `pacbio`) where there aren't many established/optimal tools yet
- **Better bioinformatic tools**: Todd's lab has a lot of tools that are fast and efficient for this. 
        - This project will be a good wrapper for promoting all the good tools produced by the lab ; and benefits of LLM
- Clear goal & **timeline**: A concrete pipeline allows us to move faster towards a paper (*outline around Sep 1st*)
 - **Competition**: Having our own custom made pipeline gives  an edge over the [highly resourced](https://seqera.io/blog/seqera-raises-26m-series-b/) seqera company's (owns nextflow) [AI-chat](https://seqera.io/ask-ai/chat) doing [similar work](https://www.healthcareittoday.com/2024/08/29/seqera-acquires-tinybio-to-advance-science-for-everyone-now-through-genai/) (news [release](https://seqera.io/blog/seqera-ai-launch/))
