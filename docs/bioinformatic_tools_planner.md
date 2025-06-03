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

## Tools list

### Quality control / preprocessing
- [NanoPlot](https://github.com/wdecoster/NanoPlot) for QC {on the raw files / final reads after filtering}
- [hostile](https://github.com/bede/hostile): To remove host contamination
- [chopper](https://github.com/wdecoster/chopper): To filter out low quality and short reads
- 


### Assembly
- Rheaa
- 

### Taxonomic classification/profiling
- Emu
- Lemur + Magnet
- Sylph
- Centrifuge

### Functional annotation
- Seqscreen
- EggnoG mapper
- Humann (? / Not for long reads - Austin? - *eukaryotic; RAM intensive*; )

### Pathogen identification
- Magnet
- Seqscreen


