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

- **Emu**: Works with example from repo. 
  
  1. Copied example from EMU repo
  2. fixed `meta` input with dummy value; solved conda issues with channel priority;
  3. Added `--db` copied from gms_16S repo
  4. Removed # -- comments from the script section and added back the blank link after `$reads` (fixed the missing input error)
  - I assume the cat .. END_VERSIONS section is for the meta outputs? This might be useful to emulate for other tools. Run gms_16S by itself with nextflow to see the output and what's going on here
    - This is present in every nf-core module. example: [fastqc](https://github.com/nf-core/modules/blob/master/modules/nf-core/fastqc/main.nf#L47)
  5. Added `--output-dir ./` to the command line arguments
    - The module was missing an output directory argument that is present in the `modules.config` of `gms_16S` repo. Adding this helps in making the module reusable without the config; other default outputs saved within results/ and not found by nextflow
    - Strange, it works fine running emu directly in the env `(/home/pbk1/micromamba/other-envs/env-ca18e434b14574137dd432ba06605711)` with `emu abundance --db databases/emu/  examples/data/emu_full_length.fa --output-dir work/emutest/` (see work/emutest)
    Error:
    ```sh
    Missing output file(s) `*abundance.tsv` expected by process `EMU_ABUNDANCE (test)`
    ```

- **Centrifuger**: Runs with the 46_1_sub10k.fastq.gz file and legionella cfr database now.
  
  - Database issue: (_figured that you need all 4 cfr_ref_idx files; hence created database from example files for the 2 legionella genomes in the example directory_)
    - Downloaded GTDB r226 index from [dropbox](https://www.dropbox.com/scl/fo/xjp5r81jxkzxest9ijxul/ADfYFKoxIyl0hrICeEI63QM?rlkey=5lij0ocrbre165pa52mavux5z&e=1&st=4ol28yv2&dl=0) ; _Downloaded only the `cfr_gtdb_r226.2.cfr` file_ ;
    - Getting same error when running directly in conda env `/home/pbk1/micromamba/other-envs/env-0ae31dcc1c5dffe0dcb45b228367f9cc` too. `centrifuger -u examples/data/46_1_sub10k.fastq.gz -x databases/cfr_gtdb_r226.2.cfr -t 4 > work/centrifugertest/test.tsv` ; _I think the database needs all 4 files_
    - Tried creating database from example files for the 2 legionella genomes in the example directory (`examples/centrifuger/ref.fa`). 

  - Still having segmentation fault error when classifying with this: Segfaults : `.command.sh: line 2: 753314 Segmentation fault      (core dumped)`
    - (ruled out) _Could be something wrong with the input file?_ Austin confirms that it works for him.. 
    - Identified that nextflow doesn't copy all the index files by just specifying the `index_prefix` path; since there are 4 separate files. _Find out how nf-core/centrifuge does it?_ : Copied their approach and it works now!
    - (Mock run works) when run directly in conda env `/home/pbk1/micromamba/other-envs/env-0ae31dcc1c5dffe0dcb45b228367f9cc` with original files works `centrifuger -u examples/data/46_1_sub10k.fastq.gz -x databases/legionella_cfr_idx/cfr_ref_idx -t 4 > work/centrifugertest/test.tsv`

  - (_Most reads are unclassified with the legionella database_) Download a mock nanopore fastq file from some nf-core module porechop etc. to test with.

## Pre-processing
- Combining : stringing pre-processing modules into a subworkflow. Testing in progress.
  - pipeline works ; but trying to incorporate default parameters that Austin wrote in as `task.ext.args` ; 
    - Switching to nf-core template for full compatibility to this one
  - (_fixed: fixed the stalling workflow with `channel.value([])` input to chopper_) Getting hung up after running runHostile now. (so the previous issue of getting stuck still persists?)
  - (_fixed: using value process `Channel.value()` for reusing single value channels_) something wrong with the runHostile subworkflow? Is only running hostile on 1 input instead of 2 ; and is getting stuck 
  - ```log
    [84/0fb496] RawNanoPlot (B011_2_sub10k)              | 2 of 2 ✔
    [b0/3ecac9] runHostile:HOSTILE_CLEAN (B011_2_sub10k) | 1 of 1 ✔
    [-        ] CHOPPER                                  -
    [-        ] FinalNanoPlot                            -
    ``` 
  - Something to do with the empty `contam_ref` channel being used up. Needs to be a value channel; [source](https://bioinformatics.stackexchange.com/questions/22161/nextflow-properly-chaining-process-outputs)
  - (_fixed: was doing single_read is 'false' instead of 'true'_) 
    - Something up with the hostile_clean running: It is in short read/bowtie2 mode (by default?) // need to add `--aligner minimap2` to the command line arguments using `task.ext.args`
    - How come it runs fine with testing then? _Might need to test the subworkflow runHostile separately?_
```log
Command error:
19:05:02 INFO: Hostile v2.0.1. Mode: paired short read (Bowtie2)
``` 
  - Not going into the `if` block (seems fixed now after adding .mmi extension)


- [x] Identify nf-core modules
- [chopper](https://nf-co.re/modules/chopper). works.
- [nanoplot](https://nf-co.re/modules/nanoplot) | [module scripts](https://github.com/nf-core/modules/blob/master/modules/nf-core/nanoplot/main.nf). Getting some issue `ModuleNotFoundError: No module named 'kaleido.scopes'`
  - test in conda env `/home/pbk1/micromamba/other-envs/env-dbb0881eb12a6a46-f15e056cde8296457e44df4488e27ca5`. Some [issue](https://github.com/wdecoster/NanoPlot/issues/417) with latest version of `python-kaleido` 1.0.
  - Tried downgrading to `python-kaleido` 0.2.1. Module runs but makes empty plots (known issue; `NanoPlot-report.html` works file)
    - Tried remaking the whole env with this dependancy locked in the `environment.yml` file (_will take care of any dependancy conflicts_)  
- [hostile_clean](https://nf-co.re/modules/hostile_clean) | [hostile_fetch](https://nf-co.re/modules/hostile_fetch). works, along with fetch (_optional_)

## nf-core compatibility
- Created a template using `nf-core pipelines create` with custom settings
  - _Assuming this is not going on nf-core since Todd would want to keep ownership rather than community owned status_
  - Still want to keep the nf-core template for composability with nf-core modules and any future forks people might make.
- Moved all components of the template `nf-core-somatem` dir into current directory and merged any similar dirs/files(`nextflow.config, modules.json, docs/, .nf-core.yml`)
  - note: extra readme saved in archive for future ideas ; config was mixed with current config for testing (_cacheDir is hardcoded path_) 


---
# Process to make a local nextflow module

1. Clone the tool's repo into a temp directory (outside this repo)
2. (optional) Test the tool with example data from the tool's own repo
3. Copy the module template from `modules/module_template.nf` to `modules/local/{tool_name}/main.nf`
4. Check for the tool's conda repo to call in the module-process's conda definition
5. For windsurf-AI's help in making the module, copy the tool's main script or readme of how to use it to the `modules/local/{tool_name}` directory as a placeholder file

Install an nf-core module using `nf-core modules install ..`


# to nf-core or not?
- modules of nf-core require a tuple input (`tuple val(meta), path(fasta)`). I was having problems with this for the emu module.

Can make the tuple using the map workflow from any pipeline_initialization subworkflow. [demo module](https://github.com/nf-core/demo/blob/1.0.2/subworkflows/local/utils_nfcore_demo_pipeline/main.nf#L75)
```groovy
.map {
            meta, fastqs ->
                return [ meta, fastqs.flatten() ]
        }
        .set { ch_samplesheet }
```

# nextflow tips
## Input files
- Need to take in files as glob patterns and create channel with metatada from them. Can use the `subworkflows/nf-core-compatibility.nf` to help with this
   - Need to use Channel.fromPath().simpleName to create meta.id from the file name
- nf-core approach seems to only take in a sample sheet and create the channel from it. If files are batched then this would be useful. 
  - Get a demo format of such an samplesheet from nf-core modules. There's the example with only id, fastq1, fastq2 columns in the default template created with `nf-core pipelines create`
- To maintain flexibility of taking in both glob patterns and sample sheet, we can copy mag's approach from [subworkflows/local/input_check.nf](https://github.com/nf-core/mag/blob/2.3.2/subworkflows/local/input_check.nf)

# Database notes

Need to record the source of each example dataset and database in the database folder here + add it to the commit message when adding any new examples? (databases won't be in the version control, maybe need a neat script that pulls them for public google drive/box.com urls)  

## Example files (`examples/`)
- `centrifuger`: Downloaded from original repo [here](https://github.com/mourisl/centrifuger/tree/master/example)
- `data/46_1_sub10k.fastq.gz` and `B01_1_sub10k.fastq.gz`: From Austin's own generated nanopore data of gut microbiome samples. Subsampled to 10k reads.
  - `data/46_1.fastq.gz`: From google drive/[example_data/agm..](https://drive.google.com/drive/u/1/folders/1MUR6sXAJSTaKXrqhVu6-rLFDw7lao5v5)
- `data/emu_full_length.fa`: From EMU repo [here](https://github.com/treangenlab/emu/tree/master/example)
- `lemur`: from original repo/[examples](https://github.com/treangenlab/lemur/tree/main/examples)
- `Sylph`: from original repo/[testfiles](https://github.com/bluenote-1577/sylph/tree/main/test_files)

## Database files (`databases/`)

### Testing/demo databases
- `Emu`: Database obtained from gms_16S repo [here](https://github.com/genomic-medicine-sweden/gms_16S/tree/master/assets/databases/emu_database)
  - Note sure if there were from the original emu? : https://osf.io/56uf7/files/osfstorage#
  - GMS-16S utilizes a combination of the ribosomal RNA Operon copy number (rrnDB) and the NCBI 16S RefSeq databases (from gms_16S [paper](https://link.springer.com/article/10.1007/s10096-025-05158-w))
- legionella_cfr_idx`: From centrifuger example files
  - mock2 test database create from example/centrifuger/ files by running `centrifuger-build -r ref.fa --taxonomy-tree nodes.dmp --name-table names.dmp --conversion-table ref_seqid.map -o ../../work/centrifugertest/legionella-cfr_ref_idx`

### Real databases

_Clean up these old notes_
- centrifuger: 
  - real database download: GTDB r226 index from [dropbox](https://www.dropbox.com/scl/fo/xjp5r81jxkzxest9ijxul/ADfYFKoxIyl0hrICeEI63QM?rlkey=5lij0ocrbre165pa52mavux5z&e=1&st=4ol28yv2&dl=0) | link derived from [centrifuger repo](https://github.com/mourisl/centrifuger#usage)
  - mock database download: [nf-core/centrifuge: minigut_cf](https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/delete_me/minigut_cf.tar.gz) | link derived from [nf-core/centrifuge](https://github.com/nf-core/modules/blob/master/modules/nf-core/centrifuge/centrifuge/tests/main.nf.test#L18C54-L18C150)
- emu: link retrieved from ?


## Updating databases, best practices
What makes certain databases automatic install from nextflow and not others?
- [mapo tofu](https://github.com/ikmb/TOFU-MAaPO)
> The pipeline can download and install the required databases for GTDBtk, MetaPhlAn and HUMAnN. Refer to the usage documentation for more details.

- [aviary](https://github.com/rhysnewell/aviary) has a nice way to determine database location (through `config`) and download them with the `--download` flag. _Check if this can use latest database or links to static versions?_

> Following tools need manual creation or download of required databases:
  - Bowtie2 (for host genome removal)
  - Kraken2 (with Braken)
  - Sylph
  - Salmon

