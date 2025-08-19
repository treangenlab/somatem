# Installation

## micromamba
Install micromamba using the command below in Linux. Source: [docs](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html)

```sh
"${SHELL}" <(curl -L micro.mamba.pm/install.sh)
```

## Nextflow
Create the micromamba environment for nextflow (and nf-core) installation : `nf_base_env`

```sh
micromamba create -f nf_base_env.yml
```

### Set the cache dir
Nextflow's cache dir is set to `~/micromamba/other-envs` in the `nextflow.config` file. So create a directory at this location the first time

```sh
mkdir -p ~/micromamba/other-envs
```

