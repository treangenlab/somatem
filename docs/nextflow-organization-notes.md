# Nextflow Organization Notes

## Rough outline of the plan
- Start with making working pipelines
- Then break the pipelines into modules
  - Use [nf-core modules](https://nf-co.re/modules) if available. Installation with `nf-core modules install <module_name>` puts them under `./modules/nf-core/modules` directory
  - If not, create a module following the [nf-core module template](https://nf-co.re/modules#module-template)
- Try to incorporate the conda env of each module/process into the nextflow definition. 
  - The path of an environment file can be specified using the conda directive [as follows](https://www.nextflow.io/docs/latest/conda.html#use-conda-environment-files):

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

## 
