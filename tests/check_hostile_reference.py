"""Run with Nextflow, nf-schema, minimap2, Hostile, Chopper and NanoPlot available."""

import csv
import gzip
import os
from pathlib import Path
import random
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[1]


def main():
    with tempfile.TemporaryDirectory(prefix="somatem-hostile-") as directory:
        tmp = Path(directory)
        rng = random.Random(94)
        host = "".join(rng.choices("ACGT", k=6000))
        reference = tmp / "host genome.fna"
        reference.write_text(f">host\n{host}\n")
        samples = ["one", "two"]
        with (tmp / "samples.csv").open("w") as sheet:
            writer = csv.writer(sheet)
            writer.writerow(["sample", "fastq_1"])
            for sample in samples:
                fastq = tmp / f"raw-{sample}.fastq.gz"
                writer.writerow([sample, fastq])
                with gzip.open(fastq, "wt") as reads:
                    for i in range(8):
                        length = 500 + i * 100
                        for name, seq in [("host", host[200:200 + length]),
                                          ("keep", "".join(rng.choices("ACGT", k=length)))]:
                            reads.write(f"@{sample}_{name}_{i}\n{seq}\n+\n{chr(65 + i) * length}\n")

        (tmp / "check.config").write_text("""
conda.enabled = false
process.resourceLimits = [cpus: 2, memory: 2.GB]
process {
    withName: '.*NanoPlot.*' { ext.args = '--no_static' }
}
""")
        subprocess.run([str(ROOT / "bin/somatem"), "pre_processing", "-i", "samples.csv",
                        "-o", "results", "-c", "check.config", "--sample_environment", "dog_gut",
                        "--hostile_reference", str(reference)], cwd=tmp, check=True,
                       env=os.environ | {"NXF_OFFLINE": "true"})
        trace_path, = (tmp / "results/pipeline_info").glob("execution_trace_*.txt")
        with trace_path.open() as trace_file:
            tasks = [row["name"] for row in csv.DictReader(trace_file, delimiter="\t")]
        assert sum("HOSTILE_INDEX" in task for task in tasks) == 1, tasks
        assert sum("HOSTILE_CLEAN" in task for task in tasks) == len(samples), tasks
        assert not any("DEACON" in task for task in tasks), tasks
        assert (tmp / "results/hostile/host genome.mmi").is_file()
        for sample in samples:
            cleaned, = [p for p in (tmp / "work").glob(f"*/*/{sample}.fastq.gz") if not p.is_symlink()]
            with gzip.open(cleaned, "rt") as reads:
                assert set(reads.read().splitlines()[::4]) == {f"@{sample}_keep_{i}" for i in range(8)}
        report = (tmp / "results/summary_reports/pre_processing_summary_report.html").read_text()
        assert "Hostile removed" in report and "reads_removed" in report
        print("PASS: two samples depleted against one host FASTA before quality filtering")


if __name__ == "__main__":
    main()
