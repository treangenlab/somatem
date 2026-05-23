# Installation

## micromamba
Install micromamba using the command below in Linux. Source: [docs](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html)

```sh
"${SHELL}" <(curl -L micro.mamba.pm/install.sh)
```

## Nextflow
Create the micromamba environment for nextflow (and nf-core) installation : `nf_base_env`. This installs the latest version of nextflow and nf-core from the `bioconda` channel.

```sh
micromamba create -f nf_base_env.yml
```

- Make sure to activate the environment before running nextflow: (each time you login)

```sh
micromamba activate nf_base_env
```

### Set the cache dir
By default, Nextflow manages Conda environments in `$workDir/conda` and creates cache directories automatically.

For power users, Somatem also supports overriding install/cache locations via environment variables before launch: `SOMATEM_HOME`, `SOMATEM_DB_DIR`, `SOMATEM_CONDA_CACHE`, and `SOMATEM_UNIFIED_DB_DIR`.


# Development notes
_This is relevant if you are modifying the pipeline/repo. Using VSCode or other IDEs based on it_

## VSCode
- Install the nextflow extension for VSCode
  - Since we are using micromamba, we need to set the nextflow > java.home path in the plugin settings (`@ext:nextflow.nextflow`) to `~/micromamba/envs/nf_base_env/lib/jvm` (within the micromamba env; **Use absolute path** by replacing `~` with your home directory)
- Install other plugins that would be useful: Rainbow csv, ?

## Cloning the repo
- Use `git clone --recurse-submodules ...` to clone the repo including it's submodules
  - If already cloned the repo the normal way, use `git submodule update --init` to update the submodules ; otherwise, the sub-module repos will be empty folders
