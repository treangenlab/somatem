"""Check both SeqScreen Bowtie2 aliases with Nextflow and a simulated killed aligner."""

import os
from pathlib import Path
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[1]


def main():
    with tempfile.TemporaryDirectory(prefix="seqscreen-bowtie2-") as directory:
        tmp = Path(directory)
        (tmp / "db").mkdir()
        (tmp / "sample.fa").write_text(">sample\nACGT\n")
        bowtie2 = tmp / "bowtie2"
        bowtie2.write_text("""#!/usr/bin/env bash
set -eu
if [[ $1 == --version ]]; then echo 'bowtie2 version test'; exit 0; fi
if [[ ${BAD_INDEX:-0} == 1 ]]; then echo 'index missing' >&2; exit 1; fi
if [[ $2 != 1 ]]; then echo '(ERR): bowtie2-align died with signal 9 (KILL)' >&2; exit 1; fi
touch "${!#}"
""")
        bowtie2.chmod(0o755)
        (tmp / "check.config").write_text(f"""
params.analysis_type = 'seqscreen'
params.threads = 32
params.max_memory_gb = 124
params.local_memory_reserve_gb = 4
includeConfig '{ROOT}/conf/base.config'
trace.enabled = true
trace.file = 'trace.txt'
trace.overwrite = true
trace.fields = 'name,attempt,exit,memory,cpus'
""")
        module = ROOT / "modules/local/seqscreen/bowtie2/main"
        (tmp / "main.nf").write_text(f"""
include {{ SEQSCREEN_BOWTIE2 as SEQSCREEN_BOWTIE2_BSAT }} from '{module}'
include {{ SEQSCREEN_BOWTIE2 as SEQSCREEN_BOWTIE2_VFDB }} from '{module}'
workflow {{
    samples = channel.of('a', 'b').map {{ id -> tuple([id: id], file('sample.fa')) }}
    SEQSCREEN_BOWTIE2_BSAT(samples, file('db'), 'bsat')
    SEQSCREEN_BOWTIE2_VFDB(samples, file('db'), 'vfdb')
}}
""")
        env = {k: v for k, v in os.environ.items() if not k.startswith("SLURM_")}
        env.update(PATH=f"{tmp}:{env['PATH']}", NXF_OFFLINE="true", BAD_INDEX="0")
        cmd = ["nextflow", "-C", "check.config", "run", "main.nf", "-ansi-log", "false"]
        subprocess.run(cmd, cwd=tmp, env=env, check=True)
        log = (tmp / ".nextflow.log").read_text()
        assert "peakRunning=1;" in log
        trace = (tmp / "trace.txt").read_text().splitlines()
        for alias in ("BSAT", "VFDB"):
            for sample in ("a", "b"):
                task = f"SEQSCREEN_BOWTIE2_{alias} ({sample}:{alias.lower()})"
                assert f"{task}\t1\t137\t124 GB\t16" in trace
                assert f"{task}\t2\t0\t124 GB\t1" in trace

        failed = subprocess.run(cmd, cwd=tmp, env=env | {"BAD_INDEX": "1"},
                                capture_output=True, text=True)
        assert failed.returncode != 0 and "index missing" in failed.stdout + failed.stderr
        assert "Execution is retried" not in (tmp / ".nextflow.log").read_text()
        print("PASS Bowtie2 local memory reservation, SIGKILL retry and ordinary failure")


if __name__ == "__main__":
    main()
