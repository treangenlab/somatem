# SeqScreen DSL2 Route

This branch adds a modular SeqScreen route selected with:

```bash
somatem seqscreen \
  --input seqscreen_samples.csv \
  --seqscreen_db /path/to/SeqScreenDB_23.4 \
  --outdir results
```

The SeqScreen samplesheet accepts either FASTA input:

```csv
sample,fasta
sample1,/path/to/sample1.fasta
```

or FASTQ input (plain or gzip-compressed):

```csv
sample,fastq
sample1,/path/to/sample1.fastq.gz
```

Input type is detected from the file extension. FASTA files continue directly into
SeqScreen. FASTQ/FASTQ.GZ files are first filtered with Chopper at a minimum mean
Phred quality of 15 (with no read-length option), then clustered with MMseqs2
`linclust` at 99% sequence identity. One representative read per cluster is
converted to FASTA and passed to the same SeqScreen workflow.

Supported modes:

```bash
--seqscreen_mode fast
--seqscreen_mode sensitive
```

The route vendors upstream SeqScreen helper scripts in `assets/seqscreen/` and runs them through DSL2 modules under `modules/local/seqscreen/`.
