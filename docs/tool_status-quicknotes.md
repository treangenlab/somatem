# Tool Status
_First test each module independently with example data from each tool's own repo_

## Taxonomic profiling
- Lemur: working with example from repo
  - Need to include the optional parameters listed in `def parse_args` function [line 79](https://github.com/treangenlab/lemur/blob/main/lemur#L79)

- Magnet: Trying to make a conda env by packaging the named dependencies in README (in the `dependencies.yml` file)
  - (_without version numbers_) python 3.9 has conflict with biopython; Need to figure out which version supports this: biopython 1.70 + conflict with other packages as well :( (experiment with fixed versions of few of them)
  - (_with version numbers_) conflicts with other named versions as well (see output below)

  ```sh
  error    libmamba Could not solve for environment specs
      The following packages are incompatible
      ├─ biopython =* * is installable with the potential options
      │  ├─ biopython [1.66|1.67|1.68|1.69|1.70] would require
      │  │  └─ python =2.7 *, which can be installed;
      │  ├─ biopython [1.66|1.67|1.68|1.69|1.70] would require
      │  │  └─ python =3.5 *, which can be installed;
      │  ├─ biopython [1.66|1.67|1.69|1.70] would require
      │  │  └─ python =3.6 *, which can be installed;
      │  ├─ biopython 1.68 would require
      │  │  └─ python =3.4 *, which can be installed;
      │  └─ biopython [1.70|1.71|...|1.85] conflicts with any installable versions previously reported;
      ├─ ete3 =3.1.2 * is not installable because it conflicts with any installable versions previously reported;
      ├─ minimap2 =2.24.r1122 * does not exist (perhaps a typo or a missing channel);
      ├─ python =3.9 * is not installable because it conflicts with any installable versions previously reported;
      └─ samtools =1.15.1 * is not installable because there are no viable options
        ├─ samtools 1.15.1 would require
        │  └─ htslib >=1.15,<1.23.0a0 * but there are no viable options
        │     ├─ htslib [1.15|1.15.1] would require
        │     │  └─ libdeflate >=1.10,<1.25.0a0 *, which conflicts with any installable versions previously reported;
        │     ├─ htslib [1.15.1|1.16|1.17] would require
        │     │  └─ libdeflate >=1.13,<1.25.0a0 *, which conflicts with any installable versions previously reported;
        │     ├─ htslib [1.17|1.18|1.19|1.19.1|1.20] would require
        │     │  └─ libdeflate >=1.18,<1.25.0a0 *, which conflicts with any installable versions previously reported;
        │     ├─ htslib [1.20|1.21] would require
        │     │  └─ libdeflate >=1.20,<1.25.0a0 *, which conflicts with any installable versions previously reported;
        │     └─ htslib [1.21|1.22] would require
        │        └─ libdeflate >=1.22,<1.25.0a0 *, which conflicts with any installable versions previously reported;
        └─ samtools 1.15.1 would require
            └─ htslib >=1.16,<1.23.0a0 *, which cannot be installed (as previously explained).
  ```

- Sylph: test module for profile with example data from repo works

- EMU: fixed `meta` input with dummy value; solved conda issues with channel priority; database not specified error
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
