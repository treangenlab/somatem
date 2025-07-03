# Tool Status
_First test each module independently with example data from each tool's own repo_

## Taxonomic profiling
- **Lemur**: working with example from repo
  - Need to include the optional parameters listed in `def parse_args` function [line 79](https://github.com/treangenlab/lemur/blob/main/lemur#L79)

- **Magnet**: conda installed. Need to make a sub-module for the dependencies in `utils` folder and test .py and .nf
  - Fixed conda env issue by channel priority (`conda-forge` before `bioconda`)

- **Sylph**: test module for profile with example data from repo works

- **EMU**: fixed `meta` input with dummy value; solved conda issues with channel priority; database not specified error
  - Copied example from EMU repo

Error:
```sh
Handling unexpected condition for
  task: name=EMU_ABUNDANCE; work-dir=null
  error [nextflow.exception.ProcessUnrecoverableException]: Not a valid path value type: groovyx.gpars.dataflow.DataflowBroadcast (DataflowBroadcast around DataflowStream[?])
```

# Process to make a nextflow module

1. Clone the tool's repo into a temp directory (outside this repo)
2. (optional) Test the tool with example data from the tool's own repo
3. Copy the module template from `modules/module_template.nf` to `modules/local/{tool_name}/main.nf`
4. Check for the tool's conda repo to call in the module-process's conda definition
5. For windsurf-AI's help in making the module, copy the tool's main script or readme of how to use it to the `modules/local/{tool_name}` directory as a placeholder file
  

# Database notes

What makes certain databases automatic install from nextflow and not others?
[mapo tofu](https://github.com/ikmb/TOFU-MAaPO)

> The pipeline can download and install the required databases for GTDBtk, MetaPhlAn and HUMAnN. Refer to the usage documentation for more details.

> Following tools need manual creation or download of required databases:
  - Bowtie2 (for host genome removal)
  - Kraken2 (with Braken)
  - Sylph
  - Salmon
