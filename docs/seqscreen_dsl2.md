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

On a shared server using the local executor, set `--threads` and `--max_memory_gb` to the CPU and RAM available to this run. For example, `--threads 8 --max_memory_gb 32` lets Nextflow schedule up to 8 CPUs and 32 GiB of task requests across all samples.

In fast mode, both DIAMOND outputs use the same block size, up to 2 billion letters. Tasks requesting less than 16 GiB use smaller blocks. Changing block size can change some alignment results. [DIAMOND memory options](https://github.com/bbuchfink/diamond/wiki/3.-Command-line-options#memory--performance-options)

The route vendors upstream SeqScreen helper scripts in `assets/seqscreen/` and runs them through DSL2 modules under `modules/local/seqscreen/`. The ONT-specific upstream path is not included in this first modular port.
