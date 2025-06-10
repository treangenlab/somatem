# Plannning tools
_This is a work in progress document to plan the tools to be included in the pipeline ([flowchart](../docs/SOMAteM-sketch-v1.2.jpg)). This should be moved to the [wiki](https://docs.github.com/en/communities/documenting-your-project-with-wikis/adding-or-editing-wiki-pages#cloning-wikis-to-your-computer) when it is finalized, for posterity. A shorter summary of the tools and what they do will be in the readme_

## Summary Table

| Category                     | Tool                                                                                                | Nextflow Implementation                                                                             |
| ---------------------------- | --------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------- |
| **Quality Control**          | [NanoPlot](https://github.com/wdecoster/NanoPlot)                                                   | [nf-core/nanoplot](https://github.com/nf-core/modules/tree/master/modules/nf-core/nanoplot)         |
|                              | [hostile](https://github.com/bede/hostile)                                                          | [nf-core/hostile](https://github.com/nf-core/modules/tree/master/modules/nf-core/hostile)           |
|                              | [chopper](https://github.com/wdecoster/chopper)                                                     | [nf-core/chopper](https://github.com/nf-core/modules/tree/master/modules/nf-core/chopper)           |
|                              | [filtlong](https://github.com/rrwick/Filtlong)                                                      | [nf-core/filtlong](https://github.com/nf-core/modules/tree/master/modules/nf-core/filtlong)         |
| **Assembly**                 | [MetaCompass](https://github.com/marbl/MetaCompass)                                                 | same: [MetaCompass](https://github.com/marbl/MetaCompass)                                           |
| *@Austin*                    | [AutoCycler](https://github.com/rrwick/AutoCycler)                                                  | -                                                                                                   |
|                              | [Flye](https://github.com/fenderglass/Flye)                                                         | [nf-core/flye](https://github.com/nf-core/modules/tree/master/modules/nf-core/flye)                 |
| **Binning**                  | ?                                                                                                   | -                                                                                                   |
| **Pangenomics**              | [tMHG-Finder](https://github.com/yongze-yin/tMHG-Finder)                                            | -                                                                                                   |
|                              | [parsnp](https://github.com/marbl/parsnp)<br>                                                       |                                                                                                     |
| **SNP/SV Detection**         | [rhea](https://github.com/treangenlab/rhea) (SV, timecourse)                                        | -                                                                                                   |
| **Taxonomic Profiling**      | [Emu](https://github.com/treangenlab/emu)                                                           | [gms-16s](https://github.com/genomic-medicine-sweden/gms_16S)                                       |
|                              | [Lemur](https://github.com/treangenlab/lemur) + [MAGnet](https://github.com/treangenlab/magnet)<br> | -                                                                                                   |
|                              | [Sylph](https://github.com/bluenote-1577/sylph)                                                     | -                                                                                                   |
|                              | [Centrifuge](https://github.com/DaehwanKimLab/centrifuge)                                           | [nf-core/centrifuge](https://github.com/nf-core/modules/tree/master/modules/nf-core/centrifuge)     |
|                              | [Centrifuger](https://github.com/mourisl/centrifuger)                                               | -                                                                                                   |
| **Functional Annotation**    | [SeqScreen](https://gitlab.com/treangenlab/seqscreen)                                               | [DSL1 to 2](https://gitlab.com/treangenlab/seqscreen/-/tree/to-dsl2?ref_type=heads)                 |
|                              | [EggNOG-mapper](https://github.com/eggnogdb/eggnog-mapper)                                          | [nf-core/eggnogmapper](https://github.com/nf-core/modules/tree/master/modules/nf-core/eggnogmapper) |
|                              | [HUMAnN](https://github.com/biobakery/humann)?                                                      | -                                                                                                   |
| **Metabolic Reconstruction** | [Bakdrive](https://gitlab.com/treangenlab/bakdrive)                                                 | -                                                                                                   |
|                              | [micom](https://github.com/micom-dev/micom)                                                         | -                                                                                                   |
|                              | [Apollo](https://genomearchitect.readthedocs.io/)                                                   | -                                                                                                   |
| **Reporting, visualization** | [FastQC](https://github.com/s-andrews/FastQC)                                                       | [nf-core/fastqc](https://github.com/nf-core/modules/tree/master/modules/nf-core/fastqc)             |
|                              | [MetagenomeScope](https://github.com/marbl/MetagenomeScope)                                         | -                                                                                                   |

Total number of unique tools: 24
Number with nextflow implmentations: 12/24

---

*Table summarizes all the tools listed below; tools are grouped by category and ordered by priority. Tools with available nextflow implementations are linked above. Tools are listed in order of priority within each category.*

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

- [NanoPlot](https://github.com/wdecoster/NanoPlot): QC plotting suite for long-read sequencing data and alignments. _for QC on the raw files / final reads after filtering_. [nf-core/module](https://github.com/nf-core/modules/tree/master/modules/nf-core/nanoplot) | [citation](https://academic.oup.com/bioinformatics/article/39/5/btad311/7160911?login=true)
- [hostile](https://github.com/bede/hostile): A tool for filtering reads that align to a host genome (_removes host contamination from microbial metagenomes_). nf-core/modules: [hostile-clean](https://nf-co.re/modules/hostile_clean/); [hostile-fetch](https://nf-co.re/modules/hostile_fetch/), [code](https://github.com/nf-core/modules/tree/master/modules/nf-core/hostile) | [citation](https://doi.org/10.1093/bioinformatics/btad728)
- [chopper](https://github.com/wdecoster/chopper): A tool to filter nanopore sequencing reads by quality and length ; *Use to filter out low quality and short reads*. [nf-core/module](https://nf-co.re/modules/chopper/), [code](https://github.com/nf-core/modules/tree/master/modules/nf-core/chopper) | [citation](https://doi.org/10.1093/bioinformatics/bty149)
- [filtlong](https://github.com/rrwick/Filtlong): Quality filtering tool for long-reads by read quality and length [nf-core/module](https://nf-co.re/modules/filtlong/), [code](https://github.com/nf-core/modules/tree/master/modules/nf-core/filtlong) | no citation


### Assembly
_@Austin working on this_

**Reference guided**
- [meta-compass](https://github.com/marbl/MetaCompass): A metagenomic reference-guided assembler that leverages multiple reference genomes

**De novo assembly**: _Note: uses a lot of RAM: 32 GB minimum_

- [Autocycler](https://github.com/rrwick/AutoCycler): Consensus long-read assembly pipeline combining multiple alternative assemblies of the same genome (using De Bruijn graph assembly) | [citation](https://www.biorxiv.org/content/10.1101/2025.05.12.653612v1.full)
- [Flye](https://github.com/fenderglass/Flye): De novo assembler for single-molecule sequencing reads using repeat graphs (PacBio and Oxford Nanopore). [nf-core/module](https://nf-co.re/modules/flye/), [code](https://github.com/nf-core/modules/tree/master/modules/nf-core/flye) | [citation](https://doi.org/10.1038/s41592-020-00971-x)
- [Canu](https://github.com/marbl/canu): Hierarchical assembler designed for high-noise single-molecule sequencing; _End of life since 2021; don't use_. [nf-core/module](https://nf-co.re/modules/canu/), [code](https://github.com/nf-core/modules/tree/master/modules/nf-core/canu) | [citation](https://doi.org/10.1101/gr.215087.116)


### Binning


**Contig based**


**Assembly graph based**

## Pangenomic analyses
- [parsnp](https://github.com/marbl/parsnp): A fast microbial core-genome alignment tool, which can output core genome phylogeny, multiple genome alignments and SNP calls. [citation](https://academic.oup.com/bioinformatics/article/40/5/btae311/7667868)  
- [tMHG-Finder](https://github.com/yongze-yin/tMHG-Finder): Tool for tree guided maximal homologous group (MHG) identification from multiple genomes. MHGs enable more accurate phylogenetic reconstruction than gene annotations, accounting for horizontal gene transfer. [citation](https://www.biorxiv.org/content/10.1101/2025.03.16.643543v1.full) & [older MHG-finder](https://doi.org/10.1371/journal.pcbi.1010216)


### SNP and SV Detection
Includes gene duplication loss

- [rhea](https://github.com/treangenlab/rhea): Detects structural variants (SV, >10 bp indels) and HGT between temporally evolving microbial metagenomic samples for large cohorts of related or similar genomes (1:n samples). [citation](https://academic.oup.com/bioinformatics/article/40/Supplement_1/i58/7700881): _Kristen et al., Bioinformatics, 2024_

### Taxonomic classification/profiling
- [Emu](https://github.com/treangenlab/emu): Taxonomic classification, and abundance estimation of 16S rRNA reads for long-read data. Nextflow DSL2 implementation: [gms-16s](https://github.com/genomic-medicine-sweden/gms_16S) + [gms-16S citation](https://link.springer.com/article/10.1007/s10096-025-05158-w) | [citation](https://www.nature.com/articles/s41592-022-01520-4)
- [Lemur](https://github.com/treangenlab/lemur): For rapid and accurate taxonomic profiling on long-read metagenomic datasets. [citation](https://www.biorxiv.org/content/10.1101/2024.06.01.596961v2.full)
	- [MAGnet](https://github.com/treangenlab/magnet): Refines taxonomic profiles for accuracy using reference genome mapping from all the reads. _same citation as Lemur_: [citation](https://www.biorxiv.org/content/10.1101/2024.06.01.596961v2.full)
- [Sylph](https://github.com/bluenote-1577/sylph): A tool for rapid and accurate species level taxonomic profiling of metagenomic data using k-mer sketches. [citation](https://www.nature.com/articles/s41587-024-02412-y): _Shaw et al., Nature Biotechnology, 2024_ | [documentation](https://sylph-docs.github.io/)
- [Centrifuge](https://github.com/DaehwanKimLab/centrifuge): A rapid and memory-efficient classification system for metagenomic sequences. [nf-core/module](https://nf-co.re/modules/centrifuge_centrifuge/), [code](https://github.com/nf-core/modules/tree/master/modules/nf-core/centrifuge) | [citation](https://doi.org/10.1101/gr.210641.116)
	- [Centrifuger](https://github.com/mourisl/centrifuger) might be better for large databases due to compression? | [citation](https://link.springer.com/article/10.1186/s13059-024-03244-4)

### Functional annotation
- [SeqScreen](https://gitlab.com/treangenlab/seqscreen): Functional screening of pathogenic sequences in metagenomic data.
[nextflow: DSL1 to 2 transition](https://gitlab.com/treangenlab/seqscreen/-/tree/to-dsl2?ref_type=heads) | [citation](https://link.springer.com/article/10.1186/s13059-022-02695-x): _Balaji et al., Genome Biology, 2022_
	- _includes **antibiotic resistance genes**_. 
- [EggNOG-mapper](https://github.com/eggnogdb/eggnog-mapper): Fast functional annotation of novel sequences using orthology assignments [nf-core/module](https://nf-co.re/modules/eggnogmapper/), [code](https://github.com/nf-core/modules/tree/master/modules/nf-core/eggnogmapper) | [citation](https://doi.org/10.1093/molbev/msab293)
- [HUMAnN](https://github.com/biobakery/humann): HMP Unified Metabolic Analysis Network - profiling microbial community metabolic potential (? / Not for long reads? (*says Austin*) - *eukaryotic; RAM intensive*; )

### Read Classification
_How is this different from taxonomic classification?_

- [SeqScreen](https://gitlab.com/treangenlab/seqscreen): Functional screening of pathogenic sequences in metagenomic data. 
- [Centrifuge](https://github.com/DaehwanKimLab/centrifuge): A rapid and memory-efficient classification system for metagenomic sequences

### Pathogen identification
- [MAGnet](https://github.com/treangenlab/magnet): Metagenomic Analysis of Genomes in the ENvironmental Toolkit
- [SeqScreen](https://gitlab.com/treangenlab/seqscreen): Functional screening of pathogenic sequences in metagenomic data

### Metabolic reconstruction
- [Bakdrive](https://gitlab.com/treangenlab/bakdrive) / [recent/private version](https://github.com/treangenlab/bakdrive): _Can take in Emu output_. [citation](https://academic.oup.com/bioinformatics/article/39/Supplement_1/i47/7210449): _Wang et al., Bioinformatics, 2023_
  - [micom](https://github.com/micom-dev/micom): _Best to use with bakdrive_. [citation](https://journals.asm.org/doi/10.1128/msystems.00606-19): _Diener et al., mSystems, 2020_
- [Apollo](https://genomearchitect.readthedocs.io/en/latest/): interactive sequence annotation editor; [citation](https://link.springer.com/article/10.1186/gb-2002-3-12-research0082). _Is this even relevant for a nextflow workflow?_

### Final: Validation / QC
_Check how MetAMOS implements this says Todd_
- Check if tools ran correctly


### Report
- [FastQC](https://github.com/s-andrews/FastQC): A quality control tool for high throughput sequence data - Can have it enhanced with LLM (_like seqera AI report does it_). [nf-core/module](https://nf-co.re/modules/fastqc/), [code](https://github.com/nf-core/modules/tree/master/modules/nf-core/fastqc) | [documentation](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
_Check how MetAMOS report was made from scratch says Todd_ 


### Visualization tools
- [MetagenomeScope](https://github.com/marbl/MetagenomeScope): Web based visualization tool for metagenomic assembly graphs. 

### Databases
- bacterial: [RefSeq](https://www.ncbi.nlm.nih.gov/refseq/)
- viral: [NCBI Viral Genomes](https://www.ncbi.nlm.nih.gov/genome/viruses/)
