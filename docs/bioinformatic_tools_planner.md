# Plannning tools
_This is a temporary document to plan the tools to be included in the pipeline ([flowchart](../docs/SOMAteM-sketch-v1.2.jpg)). This should be moved to the [wiki](https://docs.github.com/en/communities/documenting-your-project-with-wikis/adding-or-editing-wiki-pages#cloning-wikis-to-your-computer) when it is created eventually, for posterity._

Original plan from flowchart sketched on whiteboard on 2025-05-13: Inputs from Todd Treangen


## TODO:
- [ ] Add github links and citations for these tools
- [ ] List the LLM questions that will lead to each step
- [ ] Check and link the pipelines for which nextflow DSL1 or DSL2 implementation exists + citations



## Major decision points
- 3 broad branches of the pipeline
  - 16S
  - metagenomic classification
  - assembly based metagenomics
- 1 vs n samples
  - Rheaa is good for n samples so it will dictate the rest of the pipeline 
- Type of data technology: 
  - Nanopore (ONT)
  - Pacbio (Pacbio)
  - Synthetic long reads (Illumina)
- Assembly vs reads
- Reference guided vs de novo assembly

## Tools list

### Quality control / preprocessing
_@Austin working on this_

- [NanoPlot](https://github.com/wdecoster/NanoPlot) for QC {on the raw files / final reads after filtering}
- [hostile](https://github.com/bede/hostile): To remove host contamination
- [chopper](https://github.com/wdecoster/chopper): To filter out low quality and short reads
- [filtlong](https://github.com/rrwick/Filtlong): To filter out low quality and short reads


### Assembly
_@Austin working on this_

**Reference guided**
- meta-compass
- magnet?

**De novo assembly**: _Note: uses a lot of RAM: 32 GB minimum_

- Autocycler
- Flye
- Canu


### Binning
- parsnp

**Contig based**


**Assembly graph based**

## Pangenomic
- tmhg-finder


### SNP and SV Detection
Includes gene duplication loss

- Rheaa: Shines well when there's multiple samples (1:n samples)


### Metabolic reconstruction
- Bakdrive: _Can take in Emu output_
  - micom: _Best to use with bakdrive_
- apollo


### Taxonomic classification/profiling
- Emu
- Lemur + Magnet
- Sylph
- Centrifuge

### Functional annotation
- Seqscreen : Funcsocs: _includes antibiotic resistance genes_

- EggnoG mapper
- Humann (? / Not for long reads - Austin? - *eukaryotic; RAM intensive*; )

### Read Classification
_How is this different from taxonomic classification?_

- Seqscreen
- centrifuge

### Pathogen identification
- Magnet
- Seqscreen


### Final: Validation / QC
_Check how MetAMOS implements this says Todd_
- Check if tools ran correctly


### Report
- FastQC - Can have it enhanced with LLM (_like seqera AI report does it_)
_Check how MetAMOS report was made from scratch says Todd_ 


### Databases
- bacterial: RefSeq
- viral: ?  


