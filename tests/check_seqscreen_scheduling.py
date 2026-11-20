"""Run with Nextflow on PATH. Uses the SeqScreen module stubs without databases."""

import csv
import os
from pathlib import Path
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[1]


def main():
    with tempfile.TemporaryDirectory(prefix="seqscreen-scheduling-") as directory:
        tmp = Path(directory)
        (tmp / "db").mkdir()
        (tmp / "sample.fa").write_text(">sample\nACGT\n")
        (tmp / "nextflow.config").write_text(f"""
params.analysis_type = 'seqscreen'
params.threads = 32
params.max_memory_gb = 124
params.local_memory_reserve_gb = 4
params.seqscreen_evalue = 30
includeConfig '{ROOT}/conf/base.config'
process.beforeScript = {{
    def retry = task.process.endsWith('SEQSCREEN_RAPSEARCH2_BSAT') && task.attempt == 1
    // mkdir fails immediately if another stub is executing on this allocation.
    '''
    mkdir '{tmp}/running' || exit 90
    sleep 0.3
    rmdir '{tmp}/running'
    ''' + (retry ? 'exit 137' : 'true')
}}
trace.enabled = true
trace.file = 'trace.tsv'
trace.fields = 'name,attempt,exit,memory,cpus'
""")
        modules = ROOT / "modules/local/seqscreen"
        (tmp / "main.nf").write_text(f"""
include {{ SEQSCREEN_BOWTIE2 as SEQSCREEN_BOWTIE2_BSAT; SEQSCREEN_BOWTIE2 as SEQSCREEN_BOWTIE2_VFDB }} from '{modules}/bowtie2/main'
include {{ SEQSCREEN_RAPSEARCH2 as SEQSCREEN_RAPSEARCH2_BSAT; SEQSCREEN_RAPSEARCH2 as SEQSCREEN_RAPSEARCH2_VFDB }} from '{modules}/rapsearch2/main'
include {{ SEQSCREEN_CENTRIFUGE }} from '{modules}/centrifuge/main'
include {{ SEQSCREEN_DIAMOND_BLASTX }} from '{modules}/diamond_blastx/main'
workflow {{
    samples = channel.of('a', 'b').map {{ id -> tuple([id: id], file('sample.fa')) }}
    db = file('db')
    SEQSCREEN_BOWTIE2_BSAT(samples, db, 'bsat')
    SEQSCREEN_BOWTIE2_VFDB(samples, db, 'vfdb')
    SEQSCREEN_RAPSEARCH2_BSAT(samples, db, 'bsat')
    SEQSCREEN_RAPSEARCH2_VFDB(samples, db, 'vfdb')
    SEQSCREEN_CENTRIFUGE(samples, db)
    SEQSCREEN_DIAMOND_BLASTX(samples, db)
}}
""")
        env = {k: v for k, v in os.environ.items() if not k.startswith("SLURM_")}
        subprocess.run(["nextflow", "-C", "nextflow.config", "run", "main.nf", "-stub-run"],
                       cwd=tmp, env=env | {"NXF_OFFLINE": "true"}, check=True)
        with (tmp / "trace.tsv").open() as trace:
            tasks = list(csv.DictReader(trace, delimiter="\t"))
        assert len(tasks) == 14  # Six processes per sample plus two retries.
        assert all(task["memory"] == "124 GB" for task in tasks)
        assert sum(task["exit"] == "137" for task in tasks) == 2
        assert sum(task["attempt"] == "2" and task["exit"] == "0" for task in tasks) == 2
        assert all(task["cpus"] == "1" for task in tasks if task["attempt"] == "2")
        assert sum(task["exit"] == "0" for task in tasks) == 12
        print("PASS SeqScreen schedules one heavy task across tools, samples and retries")


if __name__ == "__main__":
    main()
