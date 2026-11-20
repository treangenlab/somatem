# SeqScreen DSL2 Route

This branch adds a modular SeqScreen route selected with:

```bash
somatem seqscreen \
  --input seqscreen_samples.csv \
  --seqscreen_db /path/to/SeqScreenDB_23.4 \
  --outdir results
```

The SeqScreen samplesheet uses FASTA input:

```csv
sample,fasta
sample1,/path/to/sample1.fasta
```

Supported modes:

```bash
--seqscreen_mode fast
--seqscreen_mode sensitive
```

Inside a Slurm batch job, the local executor takes its CPU budget from `--cpus-per-task` or the CPUs allocated on the current node. Somatem reads Slurm's `--mem` or `--mem-per-cpu` allocation and reserves `--local_memory_reserve_gb` (default `4`) for Nextflow and other overhead. For example, a single task with 32 CPUs and 128 GiB gives Nextflow a scheduling budget of 32 CPUs and 124 GiB when Java can see the full allocation. The command above uses these defaults automatically. You can still supply `--threads` and `--max_memory_gb` for a smaller budget. To check the budget used by a run, look for `Creating local task monitor` in `.nextflow.log`.

The local executor runs one medium- or high-resource SeqScreen task at a time across the entire batch. Bowtie2, RAPSearch2, Centrifuge, DIAMOND and the sensitive-mode searches each reserve the full detected task memory budget, including on retries. Smaller tasks such as FASTA validation can run together while the heavy tasks are idle. Searches can use multiple CPU threads. Local medium-resource searches, including Bowtie2, RAPSearch2 and Centrifuge, use one thread on retry. Somatem recognizes a signal 9 reported by the Bowtie2 and Centrifuge wrappers and a Bowtie2 abort after a C++ allocation failure (`std::bad_alloc`).

In fast mode, both DIAMOND outputs use the same block size, up to 2 billion letters. Tasks requesting less than 16 GiB use smaller blocks. Changing block size can change some alignment results. [DIAMOND memory options](https://github.com/bbuchfink/diamond/wiki/3.-Command-line-options#memory--performance-options)

The route vendors upstream SeqScreen helper scripts in `assets/seqscreen/` and runs them through DSL2 modules under `modules/local/seqscreen/`. The ONT-specific upstream path is not included in this first modular port.
