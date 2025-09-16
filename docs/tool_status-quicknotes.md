# Tool Status
_First test each module independently with example data from each tool's own repo_

## Taxonomic profiling

- Combine into a **subworkflow**: Currently have Emu for 16S and Lemur + Magnet for metagenomics with `params.data_type` key.
  - Test using `nextflow run subworkflows/local/taxonomic-profiling.nf -profile test --input_dir examples/lemur/example-data/example.fastq`
  - Note: Lemur needs full DB to run 46/B011 files ; Magnet needs > 1 hit to run clustering

- **Lemur**: working with example from repo; takes 45 m to run on 10k reads, full db (specifically for `46_1_sub10k.fastq.gz` file). (_**todo:**_ check memory requirement and give `high_memory` label? ; currently has `process_high`)
  - Tried to run `46_1_sub10k.fastq.gz` file with the full lemur database (`Refseq v221 bac..+ fungi`) and it took very long (45m, on 12 cpus, 72 GB memory). Why is the output file `abundance.tsv` so tiny? -- _is it because the reads were not cleaned?_
    ```log
    Completed at: 13-Aug-2025 17:30:31
    Duration    : 45m 46s
    CPU hours   : 9.1
    Succeeded   : 1
    ``` 
  - Made nf-core compatible (tuple input w meta, `ext.args`, `versions.yml`). (_**todo:**_ Need to expand to output files to send to MAGnet easily)
  - (_later?_) Need to include the optional parameters listed in `def parse_args` function [line 79](https://github.com/treangenlab/lemur/blob/main/lemur#L79)

- **Magnet**: Errors with ncbi datasets downloading? Debug with Eddy's Mimic project env.
  - (_debug_) Using Eddy's Mimic project env, magnet runs fine ; and there's more time gap between each entry of the downloaded genomes. **_todo_**: Look for something that is pacing the number of reqests, some other ncbi tool in conda?
    - Output log showing all 15 entries downloaded. _Stuck at the unzip step since no bash commands are found in this env_ ; Fix using `export PATH=$PATH:/usr/bin/`, this might be due to accessing an env not within the user's home directory.
    ```log
    min_abundance: 0
    num_species: 15
    751585           GCF_000210595.1         Coprococcus sp. ART55/1                 ART55/1         Chromosome
    165186           Genome Not Found.
    360807           GCF_001406855.1         Roseburia inulinivorans                 L1-83           Contig
    40520            GCF_025147765.1         Blautia obeum ATCC 29174                ATCC 291..      Complete Genome
    2981771          GCF_025567165.1         Laedolimicola ammoniilytica             Sanger_0..      Scaffold
    2292206          GCF_003478505.1         Clostridium sp. AF27-2AA                AF27-2AA        Scaffold
    1519439          GCF_000765235.1         Oscillibacter sp. ER4                   ER4             Contig
    2763663          GCF_014385265.1         Enterocloster hominis (ex Li..          BX10            Contig
    592978           GCF_051861395.1         Mediterraneibacter faecis               Bg7063          Complete Genome
    29523            Genome Not Found.
    820              GCF_044361425.1         Bacteroides uniformis                   JCM5828         Complete Genome
    853              GCF_002586945.1         Faecalibacterium prausnitzii            Indica          Complete Genome
    821              GCF_000012825.1         Phocaeicola vulgatus ATCC 84..          ATCC 848..      Complete Genome
    259315           Genome Not Found.
    745368           GCF_016900095.1         Gemmiger formicilis                     An812           Contig
    ```
      - Found a dependancy that might be relevant to the rate control of requests: `ncbi-dataset-cli` ; 
        - [ncbi-vdb](https://github.com/ncbi/ncbi-vdb) : This corresponds to SRA; not sure if this has anything to do here. 
  - Issue: `datasets download genomes ...` command not running? ; returning non-zero exit status. But it runs fine when run directly in conda env `/home/pbk1/micromamba/other-envs/env-b5d51811293ef0bc-8f384ba07576431396f94c81825f8f83` alone in bash. running through `.command.sh` fails after downloading a few genomes / not reproducible with each run..
    - When running with `check=False`, it gives `Error: 429 Too Many Requests` between a few files. 
    - error details:
    ```log
    File "/home/pbk1/micromamba/other-envs/env-b5d51811293ef0bc-8f384ba07576431396f94c81825f8f83/lib/python3.9/subprocess.py", line 528, in run
    raise CalledProcessError(retcode, process.args,
    subprocess.CalledProcessError: Command '['datasets', 'download', 'genome', 'accession', 'GCF_001406855.1', '--include', 'genome', '--filename', '46_1_sub10k.fastq-magnet-output/ncbi_downloads/360807.zip', '--no-progressbar']' returned non-zero exit status 1.
    ```
    - Progress output:
    ```log
    min_abundance: 0
    num_species: 15
    751585           GCF_000210595.1         Coprococcus sp. ART55/1                 ART55/1         Chromosome
    165186           Genome Not Found.
    360807           GCF_001406855.1         Roseburia inulinivorans                 L1-83           Contig
    40520            GCF_025147765.1         Blautia obeum ATCC 29174                ATCC 291..      Complete Genome
    2981771          GCF_025567165.1         Laedolimicola ammoniilytica             Sanger_0..      Scaffold
    ```
  - (_doesn't show up in the complex example_) Issue: `AgglomerativeClustering` , problem with `affinity` parameter; could this be because of a single thing being clustered?.
    - error details:
    ```log
    File "/home/pbk1/somatem/modules/local/magnet/magnet-repo/magnet.py", line 76, in find_representative_genome
    model = AgglomerativeClustering(affinity='precomputed', n_clusters=None, compute_full_tree=True,
    TypeError: __init__() got an unexpected keyword argument 'affinity' 
    ``` 
  - Issue: fastANI not found. (Add `fastANI` to the yml file)
  - Issue with ete3's `NCBITaxa` class: Error: `unzip: outdir/ncbi_downloads/*.zip -d outdir/ returned non-zero exit status 9.` (tested while adding `versions.yml` to the module)
    - Looks like ncbi datasets are not being downloaded hence the unzip command fails ; but when running alone, the unzip command gives message saying _no zip files have been found.._
    - added the missing dependency : `ncbi-dataset-cli` to the yml file ; still same error.
    - tried to update the ncbi taxonomy database update in python (`ncbi_taxa_db.update_taxonomy_database()`, where `ncbi_taxa_db` is an instance of `NCBITaxa` class within `ete3` package) ; This should create a +600 MB database within `~/.etetoolkit/`. The first run is supposed to do this but something might have gone wrong since I had a 47 MB file here instead.
    - taxid of `1639` which is ([listeria monocytogenes](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=1639&lvl=3&lin=f&keep=1&srchmode=5&unlock)) is still not found by `ncbi_taxa_db.get_lineage` function. 
    - (_more testing_) Testing on `46_1_sub10k.fastq.gz` file with the full lemur database classification (_has ~20 entries instead of just 1_) also fails with same error. (_seeing if more classification data fixes the issue of downloading genomes.._)
      ```sh
      nextflow run test-modules/magnet_test.nf --reads "./examples/data/46_1_sub10k.fastq.gz" --classification "./examples/lemur/46_1_sub10k.fastq-lemur-output/relative_abundance.tsv"
      ```
    - debug within dir `a1/980525cfbf580823faada78257dfd1` 
      - new env with pegged ncbi-dataset-cli at `/home/pbk1/micromamba/other-envs/env-b5d51811293ef0bc-ccf2ab982aeeb00e802d37e58f1e0dad`
      - all pagged packages: `/home/pbk1/micromamba/other-envs/env-b5d51811293ef0bc-7e8a569563b555ca18eb6d243ac1ea34`
  
  Updates:
  - conda installed. Need to make a sub-module for the dependencies in `utils` folder and test .py and .nf
  - Creating conda env for dependancies ; Fixed conda env issue by channel priority (`conda-forge` before `bioconda`)
  - Using in nextflow with a rigid [conda-lock](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#identical-conda-envs) file; build using `micromamba env export --explicit > spec-file.txt`

- **Sylph**: test module for profile with example data from repo works

- **Emu**: Works with example from repo. Copied full nf-core style from gms_16S (tuple input w meta, `ext.args`)
  
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
- Combining : stringing pre-processing modules into a subworkflow. Tested, works.
  - pipeline works ; incorporated default parameters that Austin wrote in as `task.ext.args` ; 
    - (_Switching to nf-core template for full compatibility to this one_) 
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

## Other tools
- Rhea: works with example data from repo. 
  - Couldn't handle metadata so omitted for now. _Could use the directory name as the `meta.id`?_
  - need to add outputs for each [file](https://github.com/treangenlab/rhea?tab=readme-ov-file#output-files) mentioned in the repo
  - visualization: Try [agb](https://github.com/almiheenko/AGB) for CLI visualization. _outputs to html_. Older tools include [bandage](https://github.com/rrwick/Bandage) with [cli](https://github.com/rrwick/Bandage/wiki/Command-line) option; or it's active fork [bandageNG](https://github.com/asl/BandageNG).
    - Trying agb in a separate process: `Creating env using micromamba: almiheenko::agb [cache /home/pbk1/micromamba/other-envs/env-3c441e5f6f1f6afbad3674b984213600]` ; _agb outputs a html ; I couldn't interpret the graph_
    - bandageNG: has a cli option and is on bioconda: `BandageNG image metaflye/assembly_graph.gfa bandage2.png --color meta2col.csv` ; but the 2 column csv file for colouring is different from Bandage for proper separation between contigs from different samples
    - bandage: Is on nf-core; --color options doesn't work thought in the `--help` documentation.
    - I am wondering how the intended outpyt looks like as mentioned in rhea [readme](https://github.com/treangenlab/rhea?tab=readme-ov-file#graph-visuals). 

## Orchestrating the pipeline
- Connected pre-processing, taxonomic profiling into the main workflow `somatemtem.nf`
- standardize modules: 
  1. remove separate directory from `lemur`, `magnet` modules ; output to `./` or `results/` for easier discovery? (_check other local modules too_)
  2. Make all outputs a tuple of `meta, path` except `versions.yml`
- mags plan: a) Test `somatem_mags.nf` from it's own entry workflow; 
  - [x] find example data (`input4mags`?)
  - Make databases for: `checkm2_db`, `bakta_db`, `singlem_metapackage` parameters
    ```bash
    --checkm2_db  /path/to/checkm2/uniref100.KO.1.dmnd \
    --bakta_db    /path/to/bakta/db \
    --singlem_metapackage /path/to/singlem/S5.4.0.GTDB_r226.metapackage_20250331.smpkg.zb \
    ```
  - check if `soma_test.sh` params need to be moved in or omitted (such as `--threads`, by changing `task_high` value etc.)
  - move params at the head of the workflow and simple + complex config into default location in `nextflow.config`
  - Port the entry workflow logic into the main.nf script
  - Check if the custom publishing method can be merged into the nf-core style publish? (`// Publish results to organized directories with safe copying`). Organization might be the useful case here, see how we can do it with the native publish method.
  - (_future_) : Austin will identify modules from nf-core that have been modified/moved to local eventually and add comments about changes. ([Slack](https://treangenlab.slack.com/archives/D08HP4K72QJ/p1757091054426499), 5/Sep/25) -- Could start from `somatem_mags.nf`'s diff in latest commit
    - SingleM, TaxBurst: directly in local ; Bakta moved to local ; checkm2_parse: custom made likely in local.


### Implementation notes
- Need to optimize the high memory and high threads processes : split / check the individual requirements for different processes like Austin did.
  - `Lemur`, `magnet`: 
  - `flye`: 
- Best Practices for Optimization of memory (rss) (from Seqera-AI)
  - Start Conservative: Set memory to ~1.2x the peak_rss observed
  - Monitor Efficiency: Aim for 60-80% memory utilization
  - Handle Failures: Use retry logic for memory-related failures (exit codes 137-140)
  - Consider Input Size: Scale resources based on input file sizes when possible
  ```groovy
  process SCALE_BY_INPUT {
      memory { 4.GB + (input_file.size() / 1000000000).intValue().GB }
  }
  ```
  - Use %cpu to measure efficiency of cpus - aim for 70-90% utilization

### Organization notes
_Use this opportunity of moving from `t8` to `owlet3` to make sure that the setup is fully portable and include instructions for micromamba etc. in the readme?!_


## nf-core compatibility
- Created a template using `nf-core pipelines create` with custom settings
  - _Assuming this is not going on nf-core since Todd would want to keep ownership rather than community owned status_
  - Still want to keep the nf-core template for composability with nf-core modules and any future forks people might make.
- Moved all components of the template `nf-core-somatem` dir into current directory and merged any similar dirs/files(`nextflow.config, modules.json, docs/, .nf-core.yml`)
  - note: extra readme saved in archive for future ideas ; config was mixed with current config for testing (_cacheDir is hardcoded path_) 

- `gms_16S`: pipeline is a good example of nf-core style that we can pull parts from.
  - Input as samplesheet vs dir with reads? _ex: methylseq/older pipelines takes in --input reads dir ; gms_16S takes in --input samplesheet_

- nf-core styled modules require a tuple input (`tuple val(meta), path(fasta)`). I generate such input using the helper script `subworkflows/local/utils/nf-core-compatibility.nf`. I set meta.id to the file name without the extension and meta.single_end = true for nanopore/long read data.

---
# Nextflow notes:

## Process to make a local nextflow module

1. Make a nf-core template module using `nf-core modules create`. Or for a barebones version, copy the module template from `modules/module_template.nf` to `modules/local/{tool_name}/main.nf`
2. Check for the tool's conda repo to call in the module-process's conda definition. If the tool itself doesn't exist on conda, then get all it's dependancies in the conda env and
3. Clone the tool's repo within the `modules/local/{tool_name}` directory as a submodule using `git submodule add <tool_repo_url> modules/local/{tool_name}`
4. For windsurf-AI's help in making the module, copy the tool's main script or readme of how to use it to the `modules/local/{tool_name}` directory as a placeholder file
5. Test the module with a testing workflow that gives minimal example data. Copy the template from `test-modules/` directory
6. If the module fails, try running it in the nextflow generated conda env manually with the `bash .command.sh` in the work/.. directory

If module exists on nf-core,
- Install the nf-core module using `nf-core modules install ..`


## Nextflow tips
### Input files
- Need to take in files as glob patterns and create channel with metatada from them. Can use the `subworkflows/nf-core-compatibility.nf` to help with this
   - Need to use Channel.fromPath().simpleName to create meta.id from the file name
- nf-core approach seems to only take in a sample sheet and create the channel from it. If files are batched then this would be useful. 
  - Get a demo format of such an samplesheet from nf-core modules. There's the example with only id, fastq1, fastq2 columns in the default template created with `nf-core pipelines create`
- To maintain flexibility of taking in both glob patterns and sample sheet, we can copy mag's approach from [subworkflows/local/input_check.nf](https://github.com/nf-core/mag/blob/2.3.2/subworkflows/local/input_check.nf)  


# data/databases to download
Recording the source of each example dataset and database in the database folder here + add it to the commit message when adding any new examples? (databases won't be in the version control, maybe need a neat script that pulls them for public google drive/box.com urls)  


## Example files (`examples/`)
All example files are stored in google drive/[data/examples](https://drive.google.com/drive/u/1/folders/11ZRpUCRrhdcJarlYdMSEDlCFl3oIz6Bh)
- `data/mock9_sub10k.fastq.gz`: From zymo mock data, subsampled to 10k reads using `seqtk sample -s100 /home/Users/pacbio_bakeoff/data/ZymoMockD6331/ont/SRR17913200.fastq 10000 | gzip > assets/examples/data/mock9_sub10k.fastq.gz` (_added `gzip` later_)
- `data/mock20_sub10k.fastq.gz`: From zymo mock data, subsampled to 10k reads using `seqtk sample -s100 /home/Users/pacbio_bakeoff/data/ZymoMockD6331/ont/SRR17913199.fastq 10000 | gzip > assets/examples/data/mock20_sub10k.fastq.gz`
  - Note: get original data from [SRA](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP358686&search=WGS%20AND%20GridIon&o=instrument_s%3Aa%3Bacc_s%3Aa) if needed

Other tools' example files:
- `data/emu_full_length.fa`: From EMU repo [here](https://github.com/treangenlab/emu/tree/master/example)
- `lemur`: from original repo/[examples](https://github.com/treangenlab/lemur/tree/main/examples)
- `Sylph`: from original repo/[testfiles](https://github.com/bluenote-1577/sylph/tree/main/test_files)
- `centrifuger`: Downloaded from original repo [here](https://github.com/mourisl/centrifuger/tree/master/example)
- `data/rhea`: 2 `.fasta` files from OSF.io storage/[examples](https://osf.io/fvhw8/files/osfstorage#)

### Zymo mock
Would be nice to have a [zymobiomics microbial community standards](https://www.zymoresearch.com/collections/zymobiomics-microbial-community-standards) dataset to test the pipeline with ; pick files that take a short time to run (ex: `46_1_sub10k.fastq.gz` takes 45m to run lemur; we want under 5 mins.)
- Notes: ZymoBIOMICS® Microbial Community Standard contains three easy-to-lyse bacteria, five tough-
to-lyse bacteria, and two tough-to-lyse yeasts ; [data sheet](https://files.zymoresearch.com/datasheets/ds1706_zymobiomics_microbial_community_standards_data_sheet.pdf)
  - Might be able to use reduced databases with only these 8-10 organisms (_but this will take a while to make ; so do it later_)
- Eddy has some zymo mock data here `/home/Users/pacbio_bakeoff/data/ZymoMockD6331/ont/SRR17913200.fastq` (pacbio also exists)
: This is a 54 GB file of Zymo-gut-mock-Kit9 sample ; check [details](https://trace.ncbi.nlm.nih.gov/Traces/index.html?view=run_browser&acc=SRR17913200&display=metadata) on SRA.  
  - There are other samples in this [SRA](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP358686&search=WGS%20AND%20GridIon&o=instrument_s%3Aa%3Bacc_s%3Aa) with `Library Name`s `M46 - 50` ; not sure what these mean.


### Setup automatic download script
Need a nice way to download and arrange all the example files (for testing repo). Extension: is there any benefit to making this into a nextflow process? _simplify call / add as a preinstall step_ : COuld put this in `subworkflow/local/utils/example_data_download.nf`
- Could use `curl` as suggested below or use google drive's cli tool `gdown` : [baeldung link](https://www.baeldung.com/linux/download-large-file-gdrive-cli) 
_procedure suggested by perplexity_
  - I have a script: `assets/scripts/download_gdrive.sh` that downloads files from google drive and arranges them in the correct directory structure.
- Get Google drive's direct download link from the file's shareable link in this format `https://drive.google.com/uc?export=download&id=YOUR_FILE_ID`
  - YOUR_FILE_ID is found in the Google Drive URL (e.g., in https://drive.google.com/file/d/FILE_ID/view)
  - Thoughts:
    - Use `curl` to download the file (or `wget`)
    - Use `tar` to extract the file
    - (optional) Use `mv` to move the file to the correct location
- Use a script with a text file scraped for downloading multiple files with proper naming
 - create a text/csv file in this format: `GoogleDrive,1AbcDXYZ12345,example_data1.zip`
 - bash script:
    ```bash
    #!/bin/bash

    while IFS=, read -r TYPE LINK DEST; do
      if [[ $TYPE == "GoogleDrive" ]]; then
        wget --no-check-certificate "https://drive.google.com/uc?export=download&id=${LINK}" -O "${DEST}"
      elif [[ $TYPE == "Dropbox" ]]; then
        wget --no-check-certificate "${LINK}" -O "${DEST}"
      fi
    done < files.txt
    ```
- **caveats**: Limits: For large Google Drive files (>100MB), you may need additional logic to handle Google’s virus scan/interruption page. For most small files, the above works.


## Database files (`databases/`)
Should we use shared databases from Todd's group or download our own? (For Emu, Lemur, Magnet, ..?)
Context: _Moving the repo to owlet3 for space concerns on t8's `/home`_
- Shared: benefit of space ; Store updated versions separately
  - Can make a DB_dir=`/home/dbs/` and encourage users to make a similar shared dir for the dbs and update this variable in the config file
  - Scripts: make a script to look for the required db in the DB_dir and download if not found (based on the tools being used..)
- Download: benefit of modularity ; ready to deploy on other machines ; can test the scripts easily to download the dbs.. 
- Ideas for DB scripts: [Emu: osfclient](https://github.com/treangenlab/emu?tab=readme-ov-file#1-download-database) ; 

Automatic DB download ideas:
- advanced: use the [storeDir](https://www.nextflow.io/docs/latest/reference/process.html#storedir) feature to store the db in a shared location. As mentioned in [seqera forum](https://community.seqera.io/t/prevent-nextflow-from-running-a-process-if-the-output-file-exists/1723)
- can string together a module that downloads the db and relocates it to the correct location (like runHostile subworkflow)



### Real databases
_locate or reuse databases in Todd's shared dir_ `/home/dbs/` (_to minimize redundancy_)

- hostile: using default `human-t2t-hla-argos985-mycob140` using the `hostile/fetch` module, source: [hostile readme](https://github.com/bede/hostile?tab=readme-ov-file#indexes) 
- Lemur: (dir: `/home/dbs/lemur_221_db/`) Database (RefSeq v221 bacterial and archaeal genes, and RefSeq v222 fungal genes) link mentioned in the [repo](https://github.com/treangenlab/lemur?tab=readme-ov-file#obtaining-the-database). [zenodo link](https://zenodo.org/records/10802546/files/rv221bacarc-rv222fungi.tar.gz?download=1) 
- `Emu`: Database obtained from gms_16S repo [here](https://github.com/genomic-medicine-sweden/gms_16S/tree/master/assets/databases/emu_database)
  - Note sure if there were from the original emu? : https://osf.io/56uf7/files/osfstorage#
  - GMS-16S utilizes a combination of the ribosomal RNA Operon copy number (rrnDB) and the NCBI 16S RefSeq databases (from gms_16S [paper](https://link.springer.com/article/10.1007/s10096-025-05158-w))
- checkm2_db: (dir: `/home/dbs/checkm2_db/`) : uniref100.KO.1.dmnd
  - Use the `checkm2_download` script from `nf-core/checkm2` to download the database? _the file needs to be relocated, similar to hostile fetch_
  - Make a custom script : might have issues with writing within the conda env [#51](https://github.com/chklovski/CheckM2/issues/73)

later: 
- centrifuger (_not downloaded_): GTDB r226 index from [dropbox](https://www.dropbox.com/scl/fo/xjp5r81jxkzxest9ijxul/ADfYFKoxIyl0hrICeEI63QM?rlkey=5lij0ocrbre165pa52mavux5z&e=1&st=4ol28yv2&dl=0) | link derived from [centrifuger repo](https://github.com/mourisl/centrifuger#usage)



### Testing/demo databases
- legionella_cfr_idx`: From centrifuger example files
  - mock2 test database create from example/centrifuger/ files by running `centrifuger-build -r ref.fa --taxonomy-tree nodes.dmp --name-table names.dmp --conversion-table ref_seqid.map -o ../../work/centrifugertest/legionella-cfr_ref_idx`
- centrifuger: mock database download: [nf-core/centrifuge: minigut_cf](https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/delete_me/minigut_cf.tar.gz) | link derived from [nf-core/centrifuge](https://github.com/nf-core/modules/blob/master/modules/nf-core/centrifuge/centrifuge/tests/main.nf.test#L18C54-L18C150)



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

