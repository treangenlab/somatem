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
Nextflow's cache dir is set to `~/micromamba/other-envs` in the `nextflow.config` file's `conda.cacheDir` parameter. So create a directory at this location the first time.
_The default cache dir for nextflow created environments is `Somatem/work/..`_; Since I like to delete `work/` frequently to save space, I set this to a different location outside the repo.

```sh
mkdir -p ~/micromamba/other-envs
```


# Development notes
_This is relevant if you are modifying the pipeline/repo. Using VSCode or other IDEs based on it_

## VSCode
- Install the nextflow extension for VSCode
  - Since we are using micromamba, we need to set the nextflow > java.home path in the plugin settings (`@ext:nextflow.nextflow`) to `/home/pbk1/micromamba/envs/nf_base_env/lib/jvm` (within the micromamba env)
- Install other plugins that would be useful: Rainbow csv, ?

## Cloning the repo
- Use `git clone --recurse-submodules ...` to clone the repo including it's submodules
  - If already cloned the repo the normal way, use `git submodule update --init` to update the submodules ; otherwise, the sub-module repos will be empty folders
