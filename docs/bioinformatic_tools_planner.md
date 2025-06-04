# Plannning tools
_This is a temporary document to plan the tools to be included in the pipeline ([flowchart](../docs/SOMAteM-sketch-v1.2.jpg)). This should be moved to the [wiki](https://docs.github.com/en/communities/documenting-your-project-with-wikis/adding-or-editing-wiki-pages#cloning-wikis-to-your-computer) when it is created eventually, for posterity._

Original plan from flowchart sketched on whiteboard on 2025-05-13: Inputs from Todd Treangen


## TODO:
- [x] Add github links and citations for these tools
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

- [NanoPlot](https://github.com/wdecoster/NanoPlot): QC plotting suite for long read sequencing data and alignments. _for QC on the raw files / final reads after filtering_
- [hostile](https://github.com/bede/hostile): A tool for filtering reads that align to a host genome (_removes host contamination from microbial metagenomes_)
- [chopper](https://github.com/wdecoster/chopper): A tool to filter nanopore sequencing reads by quality and length ; *Use to filter out low quality and short reads*
- [filtlong](https://github.com/rrwick/Filtlong): Quality filtering tool for long reads by removing the worst read segments


### Assembly
_@Austin working on this_

**Reference guided**
- [meta-compass](https://github.com/marbl/MetaCompass): A metagenomic reference-guided assembler that leverages multiple reference genomes
- [MAGnet](https://github.com/treangenlab/magnet): Metagenomic Analysis of Genomes in the ENvironmental Toolkit

**De novo assembly**: _Note: uses a lot of RAM: 32 GB minimum_

- [Autocycler](https://github.com/rrwick/AutoCycler): Automated long read assembly pipeline combining multiple assemblers
- [Flye](https://github.com/fenderglass/Flye): De novo assembler for single-molecule sequencing reads, such as those produced by PacBio and Oxford Nanopore
- [Canu](https://github.com/marbl/canu): A fork of the Celera Assembler designed for high-noise single-molecule sequencing (PacBio, Oxford Nanopore)


### Binning
- [parsnp](https://github.com/marbl/parsnp): A fast core-genome alignment and SNP detection tool for microbial genomes

**Contig based**


**Assembly graph based**

## Pangenomic
- [TMHG-Finder](https://github.com/treangenlab/tmhg-finder): Tool for identifying and analyzing tandem mobile genetic elements


### SNP and SV Detection
Includes gene duplication loss

- [rhea](https://github.com/treangenlab/rhea): Rapid haplotype estimation for large cohorts of related or similar genomes (1:n samples)


### Metabolic reconstruction
- [Bakdrive](https://gitlab.com/treangenlab/bakdrive) / [recent version](https://github.com/treangenlab/bakdrive): _Can take in Emu output_
  - [micom](https://github.com/micom-dev/micom): _Best to use with bakdrive_
- [Apollo](https://genomearchitect.readthedocs.io/en/latest/): interactive sequence annotation editor; [citation](https://link.springer.com/article/10.1186/gb-2002-3-12-research0082). _Is this even relevant for a nextflow workflow?_


### Taxonomic classification/profiling
- [Emu](https://github.com/treangenlab/emu): Taxonomic classification of metagenomic reads from complex microbial communities
- [Lemur](https://github.com/treangenlab/lemur): For rapid and accurate taxonomic profiling on long-read metagenomic datasets
  - [MAGNET](https://github.com/treangenlab/magnet): Metagenomic Analysis of Genomes in the ENvironmental Toolkit
- [Sylph](https://github.com/treangenlab/sylph): A tool for rapid and accurate taxonomic profiling of metagenomic data
- [Centrifuge](https://github.com/DaehwanKimLab/centrifuge): A rapid and memory-efficient classification system for metagenomic sequences

### Functional annotation
- [SeqScreen](https://gitlab.com/treangenlab/seqscreen): Functional screening of pathogenic sequences in metagenomic data
  - _includes antibiotic resistance genes_

- [EggNOG-mapper](https://github.com/eggnogdb/eggnog-mapper): Fast functional annotation of novel sequences using orthology assignments
- [HUMAnN](https://github.com/biobakery/humann): HMP Unified Metabolic Analysis Network - profiling microbial community metabolic potential (? / Not for long reads - Austin? - *eukaryotic; RAM intensive*; )

### Read Classification
_How is this different from taxonomic classification?_

- [SeqScreen](https://gitlab.com/treangenlab/seqscreen): Functional screening of pathogenic sequences in metagenomic data
- [Centrifuge](https://github.com/DaehwanKimLab/centrifuge): A rapid and memory-efficient classification system for metagenomic sequences

### Pathogen identification
- [MAGNET](https://github.com/treangenlab/magnet): Metagenomic Analysis of Genomes in the ENvironmental Toolkit
- [SeqScreen](https://gitlab.com/treangenlab/seqscreen): Functional screening of pathogenic sequences in metagenomic data


### Final: Validation / QC
_Check how MetAMOS implements this says Todd_
- Check if tools ran correctly


### Report
- [FastQC](https://github.com/s-andrews/FastQC): A quality control tool for high throughput sequence data - Can have it enhanced with LLM (_like seqera AI report does it_)
_Check how MetAMOS report was made from scratch says Todd_ 


### Databases
- bacterial: [RefSeq](https://www.ncbi.nlm.nih.gov/refseq/)
- viral: [NCBI Viral Genomes](https://www.ncbi.nlm.nih.gov/genome/viruses/)
