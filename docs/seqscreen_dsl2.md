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

The route vendors upstream SeqScreen helper scripts in `assets/seqscreen/` and runs them through DSL2 modules under `modules/local/seqscreen/`. The ONT-specific upstream path is not included in this first modular port.
