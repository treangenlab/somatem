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
| `somatem pre_processing` | NanoPlot QC, Chopper filtering, and Deacon host depletion |
| `somatem 16S` | Emu 16S profiling and Taxburst visualization |
| `somatem taxonomic_profiling` | Metagenomic profiling with Sylph (default), Lemur + MAGnet, Kraken2, or SingleM |
| `somatem assembly_mags` | Flye assembly, AGB graph reports, read mapping, Pigeon/SemiBin2 binning, CheckM2 assessment, SingleM appraisal, and Bakta annotation |
| `somatem isolate_analysis` | Autocycler/Flye isolate assembly, optional Polypolish/Pypolca hybrid polishing, Kraken2, CheckM2, MOB-suite, Bakta, and BTyper3 |
| `somatem seqscreen` | Modular SeqScreen pathogen screening; FASTQ inputs are Q15-filtered and clustered at 99% identity before analysis |
| `somatem genome_dynamics` | Rhea structural-variation and horizontal-gene-transfer analysis with Bandage output |
| `somatem pre_download_databases` | Pre-fetch databases for the selected configuration |

Each subworkflow will have a dedicated page in the [Somatem wiki](https://github.com/treangenlab/somatem/wiki) covering inputs, options, databases, and examples. Common options are `--input`, `--outdir`, `--threads`, `--data_type`, and `-params-file`; run `somatem help` for the current CLI reference.

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
