[![Anaconda-Server Badge](https://anaconda.org/bioconda/somatem/badges/version.svg)](https://anaconda.org/bioconda/somatem)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/somatem/badges/latest_release_date.svg)](https://anaconda.org/bioconda/somatem)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/somatem/badges/downloads.svg)](https://anaconda.org/bioconda/somatem)

# somatem

**A modular, open-source Nextflow toolkit for long-read microbiome analysis.**

Somatem supports Oxford Nanopore and PacBio 16S, metagenomic, isolate, and longitudinal datasets. It provides preprocessing, taxonomic profiling, assembly and MAG recovery, isolate characterization, pathogen screening, and genome-dynamics workflows with reproducible Conda environments and HTML summaries.

## Install

Somatem supports Linux and macOS and requires Nextflow 26.04.6 or newer. Install the Bioconda package with conda, mamba, or micromamba:

```bash
micromamba create -n somatem -c bioconda somatem
micromamba activate somatem
somatem help
```

Download the example datasets and sample sheets before trying a workflow without `--input`:

```bash
somatem examples --outdir /path/to/somatem-examples
somatem 16S --outdir results
```

The downloaded location is remembered automatically; set `SOMATEM_EXAMPLE_DIR` to override it.

### Environment variables

The defaults work for local installations, but shared systems and HPC users may want to relocate databases and caches:

- `SOMATEM_DB_DIR`: downloaded databases (default: `~/somatem_databases`)
- `NXF_CONDA_CACHEDIR`: generated Nextflow Conda environments
- `MAMBA_ROOT_PREFIX`: downloaded and extracted micromamba packages
- `SOMATEM_UNIFIED_DB_DIR`: unified databases used for ensemble species detection (defaults to `SOMATEM_DB_DIR`)

Set overrides in your shell or source [`assets/scripts/somatem_env.sh`](assets/scripts/somatem_env.sh) from your shell profile:

```bash
export SOMATEM_DB_DIR=/path/to/databases
echo "source /path/to/somatem_env.sh" >> ~/.bashrc
```

## Workflows

| Command | Purpose |
| --- | --- |
| `somatem pre_processing` | NanoPlot QC, Chopper filtering, and optional host removal with Deacon or Hostile |
| `somatem 16S` | Emu 16S profiling and Taxburst visualization |
| `somatem taxonomic_profiling` | Metagenomic profiling with Sylph (default), Lemur + MAGnet, Kraken2, or SingleM |
| `somatem assembly_mags` | Flye assembly, AGB graph reports, read mapping, Pigeon/SemiBin2 binning, CheckM2 assessment, SingleM appraisal, and Bakta annotation |
| `somatem isolate_analysis` | Autocycler consensus from Flye, Myloasm, metaMDBG, Canu, and Hifiasm, optional Polypolish/Pypolca hybrid polishing, Kraken2, CheckM2, MOB-suite, Bakta, and BTyper3 |
| `somatem seqscreen` | Modular SeqScreen pathogen and functional screening with TSV and portable HTML reports |
| `somatem genome_dynamics` | Rhea structural-variation and horizontal-gene-transfer analysis with Bandage output |
| `somatem pre_download_databases` | Pre-fetch databases for the selected configuration |

Each subworkflow will have a dedicated page in the [Somatem wiki](https://github.com/treangenlab/somatem/wiki) covering inputs, options, databases, and examples. Common options are `--input`, `--outdir`, `--threads`, `--data_type`, and `-params-file`; run `somatem help` for the current CLI reference.

### MAG assembly with a single host genome

Supply one host FASTA to select Hostile for all samples in the run. For example, to assemble dog gut samples:

```bash
somatem assembly_mags -i samples.csv -o results --sample_environment dog_gut --hostile_reference dog.fna
```

Somatem builds the minimap2 index once and uses Hostile to remove matching reads before Chopper and assembly. The reference can be `.fna`, `.fa`, or `.fasta`, optionally gzipped. To reuse the index in a later run, pass `--hostile_reference results/hostile/dog.mmi`. With the option omitted, Somatem uses the existing Deacon settings, which enable host removal by default for human and mouse metagenomes.

### Isolate batches and assembly recovery

If Autocycler fails or its consensus is too small, the isolate workflow passes along a candidate assembly that meets `--autocycler_min_genome_fraction` (default `0.8`) of that sample's supplied or estimated genome size. It selects the qualifying candidate closest to the expected size, breaking ties by filename. If no candidate qualifies, the sample is reported without an assembly and other samples continue. The size threshold is unchanged; fallback selection does not establish biological completeness.

Each completed selection writes `isolate_analysis/<sample>/autocycler/<sample>.autocycler.selection.tsv`, recording `consensus`, `fallback`, or `no_assembly`, the source, and the reason. The summary report shows recovery status; the length audit and Autocycler log remain alongside it. A fallback FASTA has no Autocycler consensus graph. Samples with all expected candidate assemblies proceed immediately; partial groups wait until the remaining assembler tasks finish.

Sample task failures are isolated to that sample and stage. Resource-related exits receive one retry, then the failed task is skipped; other sample task failures are skipped immediately. Check the summary and `pipeline_info/execution_trace_*.txt` for missing outputs and failed tasks, since a completed batch can contain failed samples. MOB-suite database readiness is checked once per run before parallel sample jobs start, reusing an initialized database. Shared database setup failures, invalid workflow settings, and aggregate report failures still fail the run.

When running isolate analysis or SeqScreen inside a Slurm batch job with the local executor, set the CPUs and memory once in the batch script. Somatem reads the allocation at startup, reserves 4 GiB by default, and uses the remaining budget to schedule concurrent task requests. Outside Slurm it detects the resources visible to Java. Use `--threads` or `--max_memory_gb` for a smaller budget. Actual tool memory use can exceed a request. Retain the Nextflow work directory and use `-resume` when rerunning a batch to reuse successful tasks.

## Storage and databases

Reference databases can require 100 GB or more; Bakta alone is approximately 60 GB. The first applicable run downloads missing databases for reuse.

```bash
export SOMATEM_HOME=/path/to/somatem-data
export SOMATEM_DB_DIR=/path/to/databases
export SOMATEM_CONDA_CACHE=/path/to/nextflow-conda-envs
export MAMBA_ROOT_PREFIX=/path/to/micromamba-cache
export SOMATEM_UNIFIED_DB_DIR=/path/to/unified-databases
```

Explicit variables override `SOMATEM_HOME`. The installed launcher otherwise keeps Nextflow environments and micromamba packages under `$CONDA_PREFIX/share/somatem`. Direct `nextflow run` launches use the active `$CONDA_PREFIX`, or `.somatem/conda_cache` when no Conda environment is active. On a cluster, cache directories must be accessible from every compute node.

The `assembly_mags` workflow is compute- and storage-intensive. Plan for large temporary files and substantially longer runtimes than preprocessing or profiling workflows.

## Documentation and citation

Developer notes are in [docs/dev_notes.md](docs/dev_notes.md), and software citations are collected in [CITATIONS.md](CITATIONS.md). If Somatem supports your research, please cite Somatem and the tools used by your selected workflow.

Contributions are welcome. Somatem is licensed under the [GNU GPLv3](LICENSE).
