# Testing specs 

Run tests using the config files `.._params.yml` with the preset parameters for convenience
Can use : `pilot_params.yml` for testing/iterations: to specify a path and paramers based on the sample sheets below.

Sample sheets list the paths of the files to be used for testing. These `.csv` files are included within the `assets/` directory
- `16S_params.yml` + `16S_sheet.csv` : 16S rRNA gene sequencing data for TAXONOMIC_PROFILING. _runs with `--data_type "16S" and --analysis_type "taxonomic-profiling"`_
- `mag_params.yml` + `mag_samplesheet.csv` : (_still testing..single file_) Metagenomics sequencing data (hi qual, subsample 100k reads) for ASSEMBLY_MAGS. _run with `--data_type "metagenomics" and --analysis_type "assembly"`_. Run log in issue thread 2 / #86 and google drive/test_run_logs/mag_assembly
- `mag_big_samplesheet.csv` : agm ran this ; 
- `meta_tax_params.yml` + `meta_tax_samplesheet.csv` : shallow subsampled (10k reads) metagenomics sequencing data of zymo mock for TAXONOMIC_PROFILING. _run with `--data_type "metagenomics" and --analysis_type "taxonomic-profiling"`_
- `timeseries_samples.csv` : Time series sequencing data. _run with `--data_type "metagenomics" and --analysis_type "assembly"`_

**Archived sample sheets:**.
_archiving and explaining edits that occured in commit_: 560e4024d20bce2b9bb79c3be151b95275ef239e

- `mag_big_samplesheet.csv` : (takes 8+ h) Metagenomics sequencing data (large) for ASSEMBLY_MAGS. _run with `--data_type "metagenomics" and --analysis_type "assembly"` using these files_. See issue #26 [comment here](https://github.com/treangenlab/Somatem/issues/26#issuecomment-3364217950) for more info.
```csv
sample,fastq_1
Abx,assets/examples/data/for_asm/abx_depl.fastq.gz
Veh,assets/examples/data/for_asm/veh_depl.fastq.gz
```
- `16S_sheet.csv` : 16S rRNA gene sequencing data for TAXONOMIC_PROFILING. _run with `--data_type "16S" and --analysis_type "taxonomic-profiling"`_. Previously included this one file subsampled from `SRR17913200` (10k reads): 
```csv
sample,fastq_1
zymoM95,assets/examples/data/zymoM95.fastq.gz
```

Notes:
- important mag runtime with tested datasets discussed in issues: [#26, thread 2](https://github.com/treangenlab/Somatem/issues/26#issuecomment-3364217950) and [#66](https://github.com/treangenlab/Somatem/issues/66)

---
# Modules: notes, status
_First test each module independently with example data from each tool's own repo_

## Taxonomic profiling

- Combine into a **subworkflow**: Currently have Emu for 16S and Lemur + Magnet for metagenomics with `params.data_type` key.
  - Test using `nextflow run subworkflows/local/taxonomic-profiling.nf -profile test --input_dir examples/lemur/example-data/example.fastq`
  - Note: Lemur needs full DB to run 46/B011 files ; Magnet needs > 1 hit to run clustering


### Older tools 
_Need to break these into individual `###` categories at some point for readability!_

- **Lemur**: working with example from repo; takes 45 m to run on 10k reads, full db. 
  - Tried to run old file with the full lemur database (`Refseq v221 bac..+ fungi`) and it took very long (45m, on 12 cpus, 72 GB memory). Why is the output file `abundance.tsv` so tiny? -- _is it because the reads were not cleaned?_
    ```log
    Completed at: 13-Aug-2025 17:30:31
    Duration    : 45m 46s
    CPU hours   : 9.1
    Succeeded   : 1
    ``` 
  - Made nf-core compatible (tuple input w meta, `ext.args`, `versions.yml`).
  - (_later?_) Need to include the optional parameters listed in `def parse_args` function [line 79](https://github.com/treangenlab/lemur/blob/main/lemur#L79)

- **Magnet**: Errors with ncbi datasets downloading? Debug with Eddy's Mimic project env.
  - (_update: using a later version of ncbi-dataset-cli solved the timeout issue_) Using Eddy's Mimic project env, magnet runs fine ; and there's more time gap between each entry of the downloaded genomes. Look for something that is pacing the number of reqests, some other ncbi tool in conda? / 
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

- **Emu**: Works with example from repo. Copied full nf-core style from gms_16S (tuple input w meta, `ext.args`)
  - feature integration: `taxburst`: Fails due to duplicate `Actinobacteria` for both class and phylum of Bifidobacteriales (confirmed in emu's db: `taxonomy.tsv`) ; _deleting this column makes taxburst work! : how to fix?_ ~ maybe update emu db with recent changes to phylum names?/ 
    - Not a robust solution but could run with `errorStrategy: 'ignore'` in `taxburst`? [read more](https://www.nextflow.io/docs/latest/reference/process.html#process-error-strategy)

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

### Centrifuger
**Status**: Not used in the pipeline currently. Keep as a backup tool : might be useful for the *Ensemble* analysis, with bakeoff DB from Eddy if the current 3 don't work out.  
 *Runs with the 46_1_sub10k.fastq.gz file and legionella cfr database now.*
  
  - Database issue: (_figured that you need all 4 cfr_ref_idx files; hence created database from example files for the 2 legionella genomes in the example directory_)
    - Downloaded GTDB r226 index from [dropbox](https://www.dropbox.com/scl/fo/xjp5r81jxkzxest9ijxul/ADfYFKoxIyl0hrICeEI63QM?rlkey=5lij0ocrbre165pa52mavux5z&e=1&st=4ol28yv2&dl=0) ; _Downloaded only the `cfr_gtdb_r226.2.cfr` file_ ;
    - Getting same error when running directly in conda env `/home/pbk1/micromamba/other-envs/env-0ae31dcc1c5dffe0dcb45b228367f9cc` too. `centrifuger -u examples/data/46_1_sub10k.fastq.gz -x databases/cfr_gtdb_r226.2.cfr -t 4 > work/centrifugertest/test.tsv` ; _I think the database needs all 4 files_
    - Tried creating database from example files for the 2 legionella genomes in the example directory (`examples/centrifuger/ref.fa`). 

  - Still having segmentation fault error when classifying with this: Segfaults : `.command.sh: line 2: 753314 Segmentation fault      (core dumped)`
    - (ruled out) _Could be something wrong with the input file?_ Austin confirms that it works for him.. 
    - Identified that nextflow doesn't copy all the index files by just specifying the `index_prefix` path; since there are 4 separate files. _Find out how nf-core/centrifuge does it?_ : Copied their approach and it works now!
    - (Mock run works) when run directly in conda env `/home/pbk1/micromamba/other-envs/env-0ae31dcc1c5dffe0dcb45b228367f9cc` with original files works `centrifuger -u examples/data/46_1_sub10k.fastq.gz -x databases/legionella_cfr_idx/cfr_ref_idx -t 4 > work/centrifugertest/test.tsv`

  - (_Most reads are unclassified with the legionella database_) Download a mock nanopore fastq file from some nf-core module porechop etc. to test with.

### Sylph
_test module for profile with example data from repo works_
- Made a local module for sylph_profile. Need to replace it with the official nf-core module [here](https://nf-co.re/modules/sylph_profile/)
  - (_old comment, not sure of it's relevance_) need to add an adapter module to get tax profiling output in mpa format (?). see here https://sylph-docs.github.io/sylph-tax/ 
  - [x] (_configured path but not tested yet_) Could use the GTDB database from our `/home/dbs/gtdb-r220-c200-dbv1.syldb`

**NOTE**: (6/Apr/26) Continue the sylph implementation from #73 for downloading viral, fungal and custom DBs (_using a csv file?_) in branch `sylph_i73` not merged yet. Current implementation focusing on the unified DBs and rebased part of the branch: `sylph_i73` at commit: f6c45ef 


DBs: Note from Eddy (2/Apr/26)
> Sylph specifically, running with its default db sources (but updated versions for each rather than its prebuilt versions; Sylph takes multiple dbs in case you are unaware of) may align better with its design purpose. But for people who like RefSeq, build a newer version RefSeq index for Sylph is fairly fast (about half hour)
- Is it worthwhile building a refseq DB for users: _There's the issue of periodic updates etc. right now we rely on the original tool authors to provide the up-to-date DBs_

> For viral, there's an issue I'd suggest you to look up on the repo. Just search viral and there's one talking about unexpected viral or something like that
- [ ] Find out about this viral DB issue of Sylph

### Ganon2 
- Implementing the nf-core module for [ganon_classify](https://nf-co.re/modules/ganon_classify/) using `nf-core modules install ganon/classify`
 - Do we also need? [ganon_report](https://nf-co.re/modules/ganon_report/): _What do we want as the output: some `tsv` similar to lemur?_ ; `ganon report --report-type abundance` might be the most relevant. Read [documentation](https://pirovc.github.io/ganon/reports/)
- We will get the standardized DB from Eddy's bakeoff work! (_get locally first_)
  - (*typically*) DB needs to be build locally with the latest Refseq etc. ; read [documentation](https://pirovc.github.io/ganon/default_databases/#commonly-used-sub-sets) for recommendations of commonly used subsets for DB

- Is `Ganon suitable for long reads?` asked on a GitHub issue [here](https://github.com/pirovc/ganon/issues/297)
  > ganon can be as well applied for long-reads. It's important to mind the usage of the thresholds: https://pirovc.github.io/ganon/classification/#cutoff-and-filter-rel-cutoff-rel-filter
  
  In this paper describing a similar method, ganon achieved good results with the thresholds -c 0.12 -e 0.9, which can be a good starting point. 

Info about [Eddy's unified DBs](## Bakeoff' Eddy's DBs) under # Databases files.

### Kraken2
- Implementing the nf-core module for [kraken2_kraken2](https://nf-co.re/modules/kraken2/kraken2/) using `nf-core modules install kraken2/kraken2`
- Get the standardized DB from Eddy's bakeoff work along with the parameters.
  - Figure out how to serve this DB to other somatem users : OSF ~ like Emu / (_more recent_) Zenodo ~ like lemur

## Species detection subworkflow
- Implementing the subworkflow for species detection. in progress
- Recent error, 17/Apr/26, 2:54 PM :
```sh
(nf_base_env) pbk1@owlet03:~/Somatem$ nextflow run main.nf -params-file assets/pilot_params.yml

 N E X T F L O W   ~  version 25.10.4

Launching `main.nf` [fervent_goldwasser] DSL2 - revision: 66dfc095de

Downloading databases for analysis type: species_detection
[-        ] ORCHESTRATE_SOMATEM:SOMATEM:PREPROCESSING:RawNanoPlot         -
[-        ] ORCHESTRATE_SOMATEM:SOMATEM:PREPROCESSING:CHOPPER             -
[-        ] ORCHESTRATE_SOMATEM:SOMATEM:PREPROCESSING:FinalNanoPlot       -
[-        ] ORCHESTRATE_SOMATEM:SOMATEM:SPECIES_DETECTION:SYLPH_PROFILE   -
[-        ] ORCHESTRATE_SOMATEM:SOMATEM:SPECIES_DETECTION:GANON_CLASSIFY  -
[-        ] ORCHESTRATE_SOMATEM:SOMATEM:SPECIES_DETECTION:KRAKEN2_KRAKEN2 -
Access to 'TAXONOMIC_PROFILING.out' is undefined since the workflow 'TAXONOMIC_PROFILING' has not been invoked before accessing the output attribute

 -- Check script 'workflows/somatem.nf' at line: 61 or see '.nextflow.log' file for more details

```



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
- Rhea: works with example data from repo (but is fasta files) ; tested the subworkflow with mock9 and 20 fastqs). 
  - made metadata input compatible with nf-core: Using mock metadata with "multiple" as `meta.id`
  - need to add outputs for each [file](https://github.com/treangenlab/rhea?tab=readme-ov-file#output-files) mentioned in the repo
  - visualization: (went with nf-core `bandage` for now) Try [agb](https://github.com/almiheenko/AGB) for CLI visualization. _outputs to html_. Older tools include [bandage](https://github.com/rrwick/Bandage) with [cli](https://github.com/rrwick/Bandage/wiki/Command-line) option; or it's active fork [bandageNG](https://github.com/asl/BandageNG).
    - Trying agb in a separate process: `Creating env using micromamba: almiheenko::agb [cache /home/pbk1/micromamba/other-envs/env-3c441e5f6f1f6afbad3674b984213600]` ; _agb outputs a html ; I couldn't interpret the graph_
    - bandageNG: has a cli option and is on bioconda: `BandageNG image metaflye/assembly_graph.gfa bandage2.png --color meta2col.csv` ; but the 2 column csv file for colouring is different from Bandage for proper separation between contigs from different samples
    - bandage: Is on nf-core; --color options doesn't work though in the `--help` documentation. (_not using color for now_)
    - I am wondering how the intended output looks like as mentioned in rhea [readme](https://github.com/treangenlab/rhea?tab=readme-ov-file#graph-visuals). 

# Orchestrating the pipeline
- Connected pre-processing, taxonomic profiling into the main workflow `somatemtem.nf`
- standardize modules: 
  1. remove separate directory from `lemur`, `magnet` module outputs ; output to `./` or `results/` for easier discovery? (_check other local modules too_)
  2. Make all outputs a tuple of `meta, path` except `versions.yml`
- mags plan: 
  - a) (_skipping this_) Test `somatem_mags.nf` from it's own entry workflow; 
  - [x] find example data (`input4mags`?)
  - [x] Make databases for: `checkm2_db`, `bakta_db`, `singlem_metapackage` parameters
    ```bash
    --checkm2_db  /path/to/checkm2/uniref100.KO.1.dmnd \
    --bakta_db    /path/to/bakta/db \
    --singlem_metapackage /path/to/singlem/S5.4.0.GTDB_r226.metapackage_20250331.smpkg.zb \
    ```
  - check if `soma_test.sh` params need to be moved in or omitted (such as `--threads`, by changing `task_high` value etc.)
  - move params at the head of the workflow and simple + complex config into default location in `nextflow.config`
  - [ ] Port the ~~entry~~ full mag workflow logic into the ~~main~~ `somatem_mags.nf` script
  
  - notes: how to handle the databases? Switched from an entry workflow separate run of `download_dbs` to a subworkflow in the main workflow. This enables the usage of the output channels from download_dbs in the main workflow (better than providing the static path of the databases ; for proper dependancy tracking _says seqera AI_). Will make empty channels that will be filled if a download db module is run (for switching between branches) / alternative is to mix all channels but unmixing becomes ugly. (_not using this for now_) 

  - Check if the custom publishing method can be merged into the nf-core style publish? (`// Publish results to organized directories with safe copying`). Organization might be the useful case here, see how we can do it with the native publish method.
  - (_future_) : Austin will identify modules from nf-core that have been modified/moved to local eventually and add comments about changes. ([Slack](https://treangenlab.slack.com/archives/D08HP4K72QJ/p1757091054426499), 5/Sep/25) -- Could start from `somatem_mags.nf`'s diff in latest commit
    - SingleM, TaxBurst: directly in local ; Bakta moved to local ; checkm2_parse: custom made likely in local.

- DB notes:
  - `hostile`: switching to `download_dbs` subworkflow for this. How to handle the output channel that has path only to a `referece/` directory? (_not the downloaded file explicityly_) ; Hostile clean and fetch modules might need to be recoded to be compatible with a `storeDir` based workflow?

## Implementation notes
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

**Organization notes** : _Use this opportunity of moving from `t8` to `owlet3` to make sure that the setup is fully portable and include instructions for micromamba etc. in the readme?!_


## nf-core compatibility
- Created a template using `nf-core pipelines create` with custom settings
  - _Assuming this is not going on nf-core since Todd would want to keep ownership rather than community owned status_
  - Still want to keep the nf-core template for composability with nf-core modules and any future forks people might make.
- Moved all components of the template `nf-core-somatem` dir into current directory and merged any similar dirs/files(`nextflow.config, modules.json, docs/, .nf-core.yml`)
  - note: extra readme saved in archive for future ideas ; config was mixed with current config for testing (_cacheDir is hardcoded path_) 

- `gms_16S`: pipeline is a good example of nf-core style that we can pull parts from.
  - Input as samplesheet vs dir with reads? _ex: methylseq/older pipelines takes in --input reads dir ; gms_16S takes in --input samplesheet_

- nf-core styled modules require a tuple input (`tuple val(meta), path(fasta)`). I generate such input using the helper script `subworkflows/local/utils/nf-core-compatibility.nf`. I set meta.id to the file name without the extension and meta.single_end = true for nanopore/long read data.

## Report generation

Generate html reports and standardized csv like files for feeding into omi's report generation system and LLM context for follow up chats.

- Nanoplot: _Currently using the `NanoStats.txt` file_
need to enable this option ([NanoPlot readme](https://github.com/wdecoster/NanoPlot))
```bash
--tsv_stats           Output the stats file as a properly formatted TSV.
```

# Execution 
_Typically execute the pipeline for testing at the end of day in a `screen` terminal to keep it persistant_
- Tip: Sometimes the screen session freezes and you will need to refresh it using
```bash
screen -d -r <session_id>
```


## github/somatem
- Initial execution was through `nextflow run main.nf --input_dir <dir> --data_type <type> --analysis_type <type>` 
- Later execution was through `nextflow run main.nf -params-file assets/16S_params.yml` with various params in the file. 
  - This is more reproducible (sharing just the params file) and easier to manage than command line arguments.
- For more convenience and discoverability, Todd wants to put the pipeline on bioconda: 
   - We should be able to run the pipeline with a simple command like `somatem 16S .._params.yml` : _explicit notation of what's being run rather than a vague main.nf which he doesn't like_

## bioconda/somatem

Initial steps were done by Austin with Bryce's help: 
- Need to create a bioconda recipe for the pipeline
- Need to test the recipe locally before submitting to bioconda
- Need to submit the recipe to bioconda

PK now coming back to verify that things work and make them compatible with the older github execution method as well.
Follow updates in issue: #109

Simplest test command for bioconda:
```sh
# run from ~/micromamba/envs/somatem/ dir
bin/somatem 16S -params ~/somatem/assets/16S_params.yml
```

### Notes
- Where can I find the `results/` dir for this run?
  - They are in the same place as the github run - See `mock1` and `mock2` below 
  - _is this because I ran from within this dir?_ using `(somatem) pbk1@owlet03:~/Somatem$ somatem 16S -params ~/somatem/assets/16S_params.yml`
```log
(nf_base_env) pbk1@owlet03:~/somatem/results/taxonomy$ ls -lt
total 20
drwxr-xr-x 2 pbk1 users 4096 Apr  9 20:13 mock1
drwxr-xr-x 2 pbk1 users 4096 Apr  9 20:10 mock2
drwxr-xr-x 2 pbk1 users 4096 Mar  4 13:36 mock9
drwxr-xr-x 2 pbk1 users 4096 Mar  4 13:36 mock20
drwxr-xr-x 2 pbk1 users 4096 Feb 26 19:30 zymo
```

- Both the bioconda and github versions seem to be executed by the same driver process.
  - Here's what the recent runs from `nextflow log` show :
```log
2026-04-06 14:57:51     36.2s           focused_gutenberg               OK      6c033e0449      276bf2ff-34f1-4b36-9b27-340dd37d69f2 nextflow run sylph-test.nf                                                                                                                                                                       
2026-04-06 15:11:03     23.8s           lonely_fermat                   OK      6c033e0449      96ef10d9-5239-4887-9607-cf8693e2652e nextflow run sylph-test.nf                                                                                                                                                                       
2026-04-06 15:14:02     19.1s           fabulous_pauling                OK      6c033e0449      c0c1faa9-195f-4b12-bb85-f5cb42da7f09 nextflow run sylph-test.nf                                                                                                                                                                       
2026-04-09 18:05:57     8.1s            insane_baekeland                ERR     66dfc095de      0a80ad9b-bc52-47c9-988d-461b7d6e5145 nextflow run /home/Users/pbk1/micromamba/envs/somatem/share/somatem-0.7.1/main.nf --data_type 16S --analysis_type taxonomic-profiling -params-file /home/Users/pbk1/somatem/assets/16S_params.yml
2026-04-09 19:09:12     6m 40s          special_lamport                 ERR     66dfc095de      789b5947-4b97-4775-aeb3-f4c6c5947f3c nextflow run /home/Users/pbk1/micromamba/envs/somatem/share/somatem-0.7.1/main.nf --data_type 16S --analysis_type taxonomic-profiling -params-file /home/Users/pbk1/somatem/assets/16S_params.yml
2026-04-09 20:05:02     8m 21s          disturbed_payne                 OK      66dfc095de      b71ab30b-29f8-429f-b315-2258964c7deb nextflow run /home/Users/pbk1/micromamba/envs/somatem/share/somatem-0.7.1/main.nf --data_type 16S --analysis_type taxonomic-profiling -params-file /home/Users/pbk1/somatem/assets/16S_params.yml
```


# Documentation

## metro map: nf-metro
Update: 
- Using nf-metro to make a base chart, later edited by hand in Inkscape 
versions 
- v1.3 : Shortened / less wordy, tool name in newline ; pre-processing top to bottom
- Find v1.3.5 : edited by hand in prashant's local computer / box folder!


- nf-metro: run with this ; need each branch to be labelled A -->|label| B
```bash
cd docs/somatem-docs/planning/flowcharts
micromamba activate nf-metro # activate the envelope

# render the metro map from example file
nf-metro render nf-metro/simple_pipeline.mmd -o nf-metro/simple_pipeline.svg --theme light

# render the actual metro map
nf-metro render flowchart_metro.mmd --theme light
```


_handy commands: code to add in .mmd_
```sh
# top level
%%metro title: nf-core/rnaseq
%%metro file: fastq_in | FASTQ
%%metro line: taxp | taxonomic profiling | #4CAF50

# within subgraphs
subgraph preprocessing [Pre-processing]
        %%metro exit: right | star_salmon, star_rsem, hisat2, bowtie2_salmon
        %%metro exit: bottom | pseudo_salmon, pseudo_kallisto

        %%metro entry: left | star_salmon, star_rsem, hisat2, bowtie2_salmon

        %%metro direction: TB

# between sections 
%% Inter-section edges
fastqc_filtered -->|star_salmon,star_rsem| star
...

```

- `metro validate` 
- `metro build`
- `metro publish`


# Assets / Data / Databases
**data/databases downloaded notes**
Recording the source of each example dataset and database in the database folder here + add it to the commit message when adding any new examples? (databases won't be in the version control, maybe need a neat script that pulls them for public google drive/box.com urls)  


## Example files (`examples/`)
All example files are stored in google drive/[data/examples](https://drive.google.com/drive/u/1/folders/11ZRpUCRrhdcJarlYdMSEDlCFl3oIz6Bh). `seqtk` is installed in `utils` micromamba env.

### metagenomic data 
- `data/mock9_sub10k.fastq.gz` (has <6 M reads): From zymo mock data with kit 9 (ZymoBIOMICS Gut Microbiome Standard
, 21 species across kingdoms, cat # : [D6331](https://files.zymoresearch.com/protocols/_d6331_zymobiomics_gut_microbiome_standard.pdf)), subsampled to 10k reads using `seqtk sample -s100 /home/Users/pacbio_bakeoff/data/ZymoMockD6331/ont/SRR17913200.fastq 10000 | gzip > assets/examples/data/mock9_sub10k.fastq.gz` (_added `gzip` later_)
- `data/mock20_sub10k.fastq.gz` (has <1.7 M reads) : From zymo mock data, subsampled to 10k reads using `seqtk sample -s100 /home/Users/pacbio_bakeoff/data/ZymoMockD6331/ont/SRR17913199.fastq 10000 | gzip > assets/examples/data/mock20_sub10k.fastq.gz`
  - Note: get original data from [SRA](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP358686&search=WGS%20AND%20GridIon&o=instrument_s%3Aa%3Bacc_s%3Aa) if needed. Understand what the samples mean: read paper about these samples [here](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-022-01415-8) : _Liu, Lei, et al. "Nanopore long-read-only metagenomics enables complete and high-quality genome reconstruction from mock and complex metagenomes." Microbiome 10.1 (2022): 209._
  - subsample to 50 K high quality reads. zymo `mock20_hiq50k.fastq.gz`:
  1. quality filter with chopper (for testing assembly). chopper env: `env-955b8cae971ef20a-fc555a3dbd46f0d6334849c854650578` (_Kept 1211025 reads out of 1679780 reads_)
  ```bash
  chopper --trim-approach best-read-segment --cutoff 15 -q 15 -l 500 -i /home/Users/pacbio_bakeoff/data/ZymoMockD6331/ont/SRR17913199.fastq > assets/examples/data/temp-mock20_hiq.fastq 
  ```
  2. subsample with seqtk (`utils` micromamba env)
  ```bash
  seqtk sample -s100 assets/examples/data/temp-mock20_hiq.fastq 50000 | gzip > assets/examples/data/mock20_hiq50k.fastq.gz
  ```

### 16S data
- `data/16S/mockm95_sub10k.fastq.gz` (SRR23926885) and `data/16S/mockm91_sub10k.fastq.gz` (SRR23926890): from Zymo gut microbiome, 21 community mock data D6331 from bioproject: [PRJNA804004](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP358686&search=GridIon&o=assay_type_s%3Aa%3Bacc_s%3Aa)
  - Chosen these two files as the smallest and the largest from the dataset used by Eddy for the bakeoff work. Not sure if Eddy used 16S data but the same [paper](https://link.springer.com/article/10.1186/s40168-022-01415-8) and bioproject had these under `Assay Type` = `AMPLICON`. 

To get the data for these 2 files use sra-toolkit to download and seqtk to subsample:
references: sra-tk [bioinformatics beginner tutorial](https://bioinformatics.ccr.cancer.gov/docs/bioinformatics-for-beginners-2025/Module1_Unix_Biowulf/L5/#fasterq-dump) ; 
```bash
# Activate sra-tools containing env
micromamba activate utils

# Download the data
prefetch SRR23926885 SRR23926890

# Convert to fastq
fasterq-dump --split-files SRR23926885 SRR23926890 # not sure why --split-files was suggested?

# subsample using seqtk
seqtk sample -s100 SRR23926885.fastq 10000 | gzip > mockm95_sub10k.fastq.gz
seqtk sample -s100 SRR23926890.fastq 10000 | gzip > mockm91_sub10k.fastq.gz
```

From zymo mock? 16S, subsampled to 10k reads using `seqtk sample -s100 /home/Users/pacbio_bakeoff/data/ZymoMockD6331/ont/SRR17913200.fastq 10000 | gzip > assets/examples/data/zymoM95.fastq.gz`

Other tools' example files:
- `data/emu_full_length.fa`: From EMU repo [here](https://github.com/treangenlab/emu/tree/master/example)
- `lemur`: from original repo/[examples](https://github.com/treangenlab/lemur/tree/main/examples)
- `Sylph`: from original repo/[testfiles](https://github.com/bluenote-1577/sylph/tree/main/test_files)
- `centrifuger`: Downloaded from original repo [here](https://github.com/mourisl/centrifuger/tree/master/example)
- `data/rhea`: 2 `.fasta` files from OSF.io storage/[examples](https://osf.io/fvhw8/files/osfstorage#)

Note:
Need smaller example files for faster iteration/testing of the assembly workflow. 
- Check this [tutorial](https://conmeehan.github.io/PathogenDataCourse/Worksheets/GenomeAssembly_Flye.html) for a 100 MB file from Zenodo [here](https://zenodo.org/record/4534098/files/DRR187567.fastq.bz2)
  - These are small isolate genome data for ONT available at SRA:[DRX178043](https://www.ncbi.nlm.nih.gov/sra/DRX178043) from paper [asm, 2019](https://journals.asm.org/doi/10.1128/mra.01212-19). 
  - _Would be preferable to have metagenome data or a mix with ~5 isolates?_ 
- explore larger subsamples of the zymo mock data (_known species advantage_) or `abx_depl.fastq.gz` (_many more species.._)?

### Zymo mock communities notes
Would be nice to have a [zymobiomics microbial community standards](https://www.zymoresearch.com/collections/zymobiomics-microbial-community-standards) dataset to test the pipeline with ; pick files that take a short time to run (ex: `46_1_sub10k.fastq.gz` takes 45m to run lemur; we want under 5 mins.)
- Notes: ZymoBIOMICS® Microbial Community Standard contains three easy-to-lyse bacteria, five tough-
to-lyse bacteria, and two tough-to-lyse yeasts ; [data sheet](https://files.zymoresearch.com/datasheets/ds1706_zymobiomics_microbial_community_standards_data_sheet.pdf)
  - Might be able to use reduced databases with only these 8-10 organisms (_but this will take a while to make ; so do it later_)
- Eddy has some zymo mock data here `/home/Users/pacbio_bakeoff/data/ZymoMockD6331/ont/SRR17913200.fastq` (pacbio also exists)
: This is a 54 GB file of Zymo-gut-mock-Kit9 sample ; check [details](https://trace.ncbi.nlm.nih.gov/Traces/index.html?view=run_browser&acc=SRR17913200&display=metadata) on SRA.  
  - There are other samples in this [SRA](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP358686&search=WGS%20AND%20GridIon&o=instrument_s%3Aa%3Bacc_s%3Aa) with `Library Name`s `M46 - 50` ; not sure what these mean.
  - This data is from a more complex mock community: (ZymoBIOMICS Gut Microbiome Standard
, 21 species across kingdoms, cat # : [D6331](https://files.zymoresearch.com/protocols/_d6331_zymobiomics_gut_microbiome_standard.pdf))

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

- [ ] Idea: 25/Feb/26 : Let's use box.com instead of google since we can use rclone to upload more easily into it! Then need to look for alternative to gdown that works for users without login/authorization (_making it public directory.._)

## Database files (`databases/`)

Older notes: _Question is how do we handle the databases in the config file of the pipeline repo?_
Should we use shared databases from Todd's group or download our own? (For Emu, Lemur, Magnet, ..?)
Context: _Moving the repo to owlet3 for space concerns on t8's `/home`_
- Shared: benefit of space ; Store updated versions separately
  - Can make a DB_dir=`/home/dbs/` and encourage users to make a similar shared dir for the dbs and update this variable in the config file
  - Scripts: make a script to look for the required db in the DB_dir and download if not found (based on the tools being used..)
- Download: benefit of modularity ; ready to deploy on other machines ; can test the scripts easily to download the dbs.. 
- Ideas for DB scripts: [Emu: osfclient](https://github.com/treangenlab/emu?tab=readme-ov-file#1-download-database) ; 

(_already implemented now_) Automatic DB download ideas:
- advanced: use the [storeDir](https://www.nextflow.io/docs/latest/reference/process.html#storedir) feature to store the db in a shared location. As mentioned in [seqera forum](https://community.seqera.io/t/prevent-nextflow-from-running-a-process-if-the-output-file-exists/1723)
- can string together a module that downloads the db and relocates it to the correct location (like runHostile subworkflow) or directly cd into the dir in the sh script (like MetaPhlAn in [mapo tofu](https://github.com/ikmb/TOFU-MAaPO))


### DB's plan

**Consistent config file for db paths:**
- Let's make a db_paths config file that stores the paths to the dbs for each tool
- There will be two versions of this file: 
  - `db_paths_default.yaml`: makes dummy paths for each tool's db
  - `db_paths_treangen.yaml`: stores the paths to the dbs for each tool in the shared dir (`/home/dbs/`)

**Automatic DB download:**
- Note: If the downloader process outputs files but the using process requires a folder as input, you can use a staging process that inputs files, outputs a folder ; this channel can feed the using process. Look for examples in the `lemur` and `emu` modules .  


### Real databases
_locate or reuse databases in Todd's shared dir_ `/home/dbs/` (_to minimize redundancy_)

- hostile: using default `human-t2t-hla-argos985-mycob140` using the `hostile/fetch` module, source: [hostile readme](https://github.com/bede/hostile?tab=readme-ov-file#indexes) 
- Lemur: (dir: `/home/dbs/lemur_221_db/`) Database (RefSeq v221 bacterial and archaeal genes, and RefSeq v222 fungal genes) link mentioned in the [repo](https://github.com/treangenlab/lemur?tab=readme-ov-file#obtaining-the-database). [zenodo link](https://zenodo.org/records/10802546/files/rv221bacarc-rv222fungi.tar.gz?download=1) 
- `Emu`: using `/home/dbs/emu/emu_db2023` dir from Eddy (49,243 sequences in species_taxid.fasta); _this is smaller than the original emu db (`emu_db_copy_from_kdc10`) mentioned in the paper from 2021 (49,301 sequences) stored in: https://osf.io/56uf7/overview (osfstorage/emu-prebuilt/emu.tar)_
  - (old) Database obtained from gms_16S repo [here](https://github.com/genomic-medicine-sweden/gms_16S/tree/master/assets/databases/emu_database)
  - note: GMS-16S utilizes a combination of the ribosomal RNA Operon copy number (rrnDB) and the NCBI 16S RefSeq databases (from gms_16S [paper](https://link.springer.com/article/10.1007/s10096-025-05158-w))
- checkm2_db: (dir: `/home/dbs/checkm2_db/`) : uniref100.KO.1.dmnd. Downloaded using `subworkflows/local/download_dbs.nf` from [zenodo](https://zenodo.org/records/14897628)
- bakta_db: (dir: `/home/dbs/bakta_db/`) : Downloaded using `subworkflows/local/download_dbs.nf` from [zenodo](https://zenodo.org/records/14916843)
- singlem_db: (dir: `/home/dbs/singlem_db/`) : Downloaded using `subworkflows/local/download_dbs.nf` from [zenodo](https://zenodo.org/records/15232972)
- sylph_gtdb: .. 

### Ensemble: Eddy's unified DBs from Bakeoff
_Use local DBs from the location directly instead of copying them_
Will implement a download module to fetch these databases from the remote location when they are uploaded to Zenodo etc. by Eddy eventually before her publication

(_abandoned: use local to save space_) Downloading the 3 key DBs using helper script `archive/scripts/copy_unified_dbs.sh`
Copies these 3 directories from `/home/Users/pacbio_bakeoff/data/ref_db/refseq03032025/` to `assets/databases/`:
- sylph_abf_030325 : 8.5 GB
- ganon2_abvf_030325 : 11 GB
- k2_abfv_030325 : +106 GB 

For more details: look for Eddy's unified databases (DB) here: 
```bash
cd /home/Users/pacbio_bakeoff/data/ref_db/refseq03032025/ # sylph_abf_030325
```

Information about each of these tools is included in the `/home/Users/pacbio_bakeoff/doc/README` file. Excerpt pasted below
```bash
windsurf /home/Users/pacbio_bakeoff/doc/README # or any text viewer..
```
_from the above document:_
- Database comparison: Prokaryotic(Arachea, Bacteria), Virus, Fungi, Human
    - Kraken2: 
        default: 255363 entries in seqid2taxid.map
        unified: 148260 entries in seqid2taxid.map
    - Centrifuge:
        default: 31.3G hpvc refseq + covid2 from genbank
        unified: 128G pvf 
    - Centrifuger:
        default: 42G hpvc refseq + covid2 from genbank
        unified: 68G abv
    - Sourmash:
        default: 79.43G abvf, 1214995 entries in lineage file
        unified: 1.8G abvf, 69095 entries in lineage file
    - Sylph:
        default: 14G gtdb-r220-c200, 113,104 species representative genomes
        unified: 8.5G abvf

- Information on building this seems to be included in this other README: `/home/Users/pacbio_bakeoff/data/ref_db/README`



### Testing/demo databases
- legionella_cfr_idx`: From centrifuger example files
  - mock2 test database create from example/centrifuger/ files by running `centrifuger-build -r ref.fa --taxonomy-tree nodes.dmp --name-table names.dmp --conversion-table ref_seqid.map -o ../../work/centrifugertest/legionella-cfr_ref_idx`
- centrifuger: mock database download: [nf-core/centrifuge: minigut_cf](https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/delete_me/minigut_cf.tar.gz) | link derived from [nf-core/centrifuge](https://github.com/nf-core/modules/blob/master/modules/nf-core/centrifuge/centrifuge/tests/main.nf.test#L18C54-L18C150)


### archive: future DBs?
- centrifuger (_not downloaded_): GTDB r226 index from [dropbox](https://www.dropbox.com/scl/fo/xjp5r81jxkzxest9ijxul/ADfYFKoxIyl0hrICeEI63QM?rlkey=5lij0ocrbre165pa52mavux5z&e=1&st=4ol28yv2&dl=0) | link derived from [centrifuger repo](https://github.com/mourisl/centrifuger#usage)

_notes from Austin: Aug 9th 2025_
- SingleM data download [url](https://wwood.github.io/singlem/tools/data) can be downloaded by running the `singlem data --output-directory /path/to/dbs/singlem`
- checkm2 database can be downloaded by running `checkm2 database --download --path /custom/path/`
- bakta database can be downloaded by running `bakta_db download --output <output-path> --type [light|full]` (this is the best method) but you can download from their zenodo archive.

### Updating databases, best practices
What makes certain databases automatic install from nextflow and not others?
- [mapo tofu](https://github.com/ikmb/TOFU-MAaPO)
> The pipeline can download and install the required databases for GTDBtk, MetaPhlAn and HUMAnN. Refer to the [db management](https://github.com/ikmb/TOFU-MAaPO?tab=readme-ov-file#database-management) section for more details.
  - example: MetaPhlAn: uses a `--updatemetaphlan` flag to update the database or download it first time using [module](https://github.com/ikmb/TOFU-MAaPO/blob/master/modules/metaphlan.nf) with `wget` after `cd ${params.metaphlan_db}` directly in the desired directory. 

- [aviary](https://github.com/rhysnewell/aviary) has a nice way to determine database location (through `config`) and download them with the `--download` flag. _Check if this can use latest database or links to static versions?_

> Following tools need manual creation or download of required databases:
  - Bowtie2 (for host genome removal)
  - Kraken2 (with Braken)
  - Sylph
  - Salmon



---
_Notes for development / maintenance_
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

