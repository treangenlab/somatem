# Nextflow Organization Notes

## Rough outline of the plan
- Start with making working pipelines as modules. Test module in isolation with example data
- Then connect the module to previous modules in a workflow or sub-workflow
  - Use [nf-core modules](https://nf-co.re/modules) if available. Installation with `nf-core modules install <module_name>` puts them under `./modules/nf-core/modules` directory
  - If not, create a module following the [nf-core module template](https://nf-co.re/modules#module-template)
- Try to incorporate the conda env of each module/process into the nextflow definition. 
  - Will start with a conda env with a pegged version of the package from bioconda etc.
  - Later we can use conda-lock files or apptainer containers for further reproducibility
  - side note: The path of an environment file can be specified using the conda directive [as follows](https://www.nextflow.io/docs/latest/conda.html#use-conda-environment-files):

### Conda vs apptainer
- conda in nextflow: *[enable](https://www.nextflow.io/docs/latest/conda.html#enabling-conda-environment) with nextflow.config directive: `conda.enabled = True`*. Extra: use [conda lockfiles](https://www.nextflow.io/docs/latest/conda.html#conda-lock-files) for faster env without resolution. 
	- can use [mamba](https://www.nextflow.io/docs/latest/conda.html#use-mamba-to-resolve-packages) for faster env resolves ; can use [existing env](https://www.nextflow.io/docs/latest/conda.html#use-existing-conda-environments) ~ no too reproducible ;  
	- Best Practices for Your Workflow (*seqera AI*)
	1. **Development phase**: Use mamba with environment.yml files for flexibility
	2. **Production/publication**: Switch to conda lock files for reproducibility
	3. **Use Wave**: Enable Wave to automatically generate both containers and lock files
	4. **Version pinning**: Pin critical tool versions in your environment.yml files


### How do we trade-off functioning vs updating over time?
Given that people will use this without contined supervision/development we are looking at tradeoff between continuous release vs periodic updates?

- We can use conda-lock files/apptainer containers to lock the versions of the packages but that prevents things from updating automatically
 - will need some manual update or some update script to go through and update the conda-lock files/apptainer containers. If this can be triggered by LLM that would be nice!
 - Same principle should apply to databases as well



```nextflow
process hello {
  conda '/some/path/my-env.yaml'

  script:
  """
  your_command --here
  """
}
```

  - The idea is that we will test all the conda envs but an end user/in a new environment will only need to have the base env installed and the nextflow will install the rest of the conda envs as needed when the process is run



## Pipeline Organization
- 



### Configuration
_Here's an response generated from Seqera AI for nextflow.config vs yaml file choice_


Based on the Nextflow documentation, here's the recommended approach for parameter configuration in your Nextflow pipeline:
Use nextflow.config for most parameters - it's the standard approach

The best practice is to use nextflow.config for your pipeline parameters. Here's why and how:
Primary Recommendation: nextflow.config

```nextflow
// nextflow.config
params {
    input = "data/*.fastq.gz"
    outdir = "results"
    genome = "GRCh38"
    max_memory = "128.GB"
    max_cpus = 16
}

process {
    executor = 'local'
    cpus = 2
    memory = '4.GB'
}
```

YAML/JSON Parameter Files - For Runtime Flexibility

However, Nextflow also supports YAML and JSON parameter files using the -params-file option, which is excellent for:

    Different experimental conditions
    Multiple dataset configurations
    Sharing parameter sets with collaborators

```yaml
params.yaml:
input: "data/experiment1/*.fastq.gz"
outdir: "results/experiment1"
genome: "GRCh38"
skip_qc: false
```

Usage:
```sh
nextflow run main.nf -params-file params.yaml
```

**Recommended Hybrid Approach**

The most flexible approach combines both:

    Default parameters in nextflow.config - for standard pipeline settings
    YAML parameter files - for experiment-specific overrides

```nextflow
// nextflow.config - defaults
params {
    // Default values
    input = null
    outdir = "results"
    genome = "GRCh38"
    skip_qc = false
    max_memory = "128.GB"
    max_cpus = 16
}

// Process configuration
process {
    executor = 'local'
    cpus = 2
    memory = '4.GB'
}
```

Then use YAML files for specific runs:
```yaml
# experiment1.yaml
input: "data/case_samples/*.fastq.gz"
outdir: "results/case_study"
skip_qc: true
```

Configuration Priority Order

Nextflow applies configuration in this priority order (highest to lowest):

    Command line parameters (--param value)
    Parameters from -params-file
    Config files specified with -c
    nextflow.config in current directory
    nextflow.config in pipeline directory
    Default values in the script

Best Practices Summary

    Use nextflow.config for pipeline defaults and process configuration
    Use YAML/JSON files for experiment-specific parameters
    Keep sensitive or frequently changing parameters in separate parameter files
    Document your parameters clearly in both formats

This approach gives you the best of both worlds: standard Nextflow configuration practices with the flexibility of external parameter files for different experimental conditions.
