# Tool Status
_First test each module independently with example data from each tool's own repo_

## Taxonomic profiling
- **Lemur**: working with example from repo
  - Need to include the optional parameters listed in `def parse_args` function [line 79](https://github.com/treangenlab/lemur/blob/main/lemur#L79)

- **Magnet**: Almost works ; need to test with a proper example with WGS fastq data?
  - conda installed. Need to make a sub-module for the dependencies in `utils` folder and test .py and .nf
  - Creating conda env for dependancies ; Fixed conda env issue by channel priority (`conda-forge` before `bioconda`)
  - Using in nextflow with a rigid [conda-lock](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#identical-conda-envs) file; build using `micromamba env export --explicit > spec-file.txt`

- **Sylph**: test module for profile with example data from repo works

- **Emu**: Missing query/-d/something error. (5/Jul/25). Can't find this option in the [EMU repo](https://github.com/treangenlab/emu?tab=readme-ov-file#abundance-estimation-parameters). 
  - This is strange:: Works fine running emu directly in the env `(/home/pbk1/micromamba/other-envs/env-ca18e434b14574137dd432ba06605711)` with `emu abundance --db databases/emu/  examples/data/emu_full_length.fa --output-dir work/emutest/` (see work/emutest)
  1. Copied example from EMU repo
  2. fixed `meta` input with dummy value; solved conda issues with channel priority;
  3. Added `--db` copied from gms_16S repo
  4. Removed # -- comments from the script section and added back the blank link after `$reads` (fixed the missing input error)
  - I assume the cat .. END_VERSIONS section is for the meta outputs? This might be useful to emulate for other tools. Run gms_16S by itself with nextflow to see the output and what's going on here

Error:
```sh
Missing output file(s) `*abundance.tsv` expected by process `EMU_ABUNDANCE (test)`
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
