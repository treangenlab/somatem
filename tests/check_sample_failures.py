#!/usr/bin/env python3
"""Run with `python3 tests/check_sample_failures.py`; requires Nextflow and Java.

Exercise the real base config using tiny tasks, without scientific tools or data.
"""

import csv
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1]


def main():
    nextflow = shutil.which("nextflow")
    assert nextflow, "Nextflow and Java must be on PATH"
    with tempfile.TemporaryDirectory(prefix="somatem-failure-check-") as tmp:
        tmp = Path(tmp)
        (tmp / "bin").mkdir()
        (tmp / "bin/mob_init").write_text(f'#!{sys.executable}\n' + '''
import os
from pathlib import Path
import sys
db = Path(sys.argv[sys.argv.index('--database_directory') + 1])
with (db / 'init_calls.txt').open('a') as handle:
    handle.write('init\\n')
(db / 'clusters.txt').write_text('mock clusters\\n')
if os.environ.get('SOMATEM_MOB_INCOMPLETE') != '1':
    (db / 'status.txt').write_text('complete\\n')
''')
        (tmp / "bin/mob_recon").write_text(f'#!{sys.executable}\n' + '''
from pathlib import Path
import sys
if '--version' in sys.argv:
    print('mob_suite test')
else:
    db = Path(sys.argv[sys.argv.index('--database_directory') + 1])
    assert (db / 'status.txt').read_text().strip() == 'complete'
    output = Path(sys.argv[sys.argv.index('--outdir') + 1])
    output.mkdir()
    (output / 'contig_report.txt').write_text('mock report\\n')
''')
        for tool in (tmp / "bin").iterdir():
            tool.chmod(0o755)
        (tmp / "task.nf").write_text('''
process SAMPLE_TASK {
    tag "$sample"
    input:
    tuple val(sample), val(exit_code)
    output:
    tuple val(sample), path("attempt.txt")
    script:
    def status = sample == 'transient' && task.attempt > 1 ? 0 : exit_code
    """
    echo ${task.attempt} > attempt.txt
    exit ${status}
    """
}
''')
        (tmp / "check.nf").write_text(
            f"include {{ MOB_SUITE_DATABASE_INIT; MOB_SUITE_RECON }} from '{ROOT}/modules/local/mob_suite/recon/main'\n" + '''
include { SAMPLE_TASK } from './task'
include { SAMPLE_TASK as ISOLATE_ANALYSIS_SUMMARY_REPORT } from './task'
include { SAMPLE_TASK as DATABASEDOWNLOAD } from './task'

workflow ISOLATE_ANALYSIS {
    if (params.scenario.startsWith('mob_')) {
        MOB_SUITE_DATABASE_INIT(file(params.mob_db).toString())
        MOB_SUITE_RECON(channel.of('alpha', 'beta', 'gamma', 'delta')
            .map { sample -> tuple([id: sample], file(params.assembly)) }, MOB_SUITE_DATABASE_INIT.out.db)
    } else if (params.scenario == 'summary') {
        ISOLATE_ANALYSIS_SUMMARY_REPORT(channel.of(tuple('summary', 1)))
    } else {
        SAMPLE_TASK(channel.of(tuple('healthy', 0), tuple('permanent', 1),
                               tuple('transient', 137), tuple('exhausted', 137)))
        SAMPLE_TASK.out.map { sample, attempt -> "$sample:${attempt.text.trim()}" }
            .collectFile(name: 'completed.txt', storeDir: params.outdir, newLine: true)
    }
}

workflow OTHER_ANALYSIS {
    SAMPLE_TASK(channel.of(tuple('other', 1)))
}

workflow SOMATEM {
    if (params.scenario == 'database') {
        DATABASEDOWNLOAD(channel.of(tuple('database', 1)))
    } else if (params.scenario == 'other') {
        OTHER_ANALYSIS()
    } else {
        ISOLATE_ANALYSIS()
    }
}

workflow { SOMATEM() }
''')
        config = tmp / "check.config"
        config.write_text(f'''
params.analysis_type = 'isolate-analysis'
params.threads = 2
params.max_cpus = null
params.max_memory_gb = 1
params.local_memory_reserve_gb = 4
includeConfig '{ROOT}/conf/base.config'
process.cpus = 1
process.memory = '128 MB'
process.executor = 'local'
process {{
    withName:MOB_SUITE_DATABASE_INIT {{ cpus = 1; memory = '128 MB' }}
    withName:MOB_SUITE_RECON {{ cpus = 1; memory = '128 MB' }}
}}
conda.enabled = false
params.publish_dir_mode = 'copy'
trace.enabled = true
trace.overwrite = true
trace.file = 'trace.txt'
trace.fields = 'process,tag,status,exit,attempt'
''')
        (tmp / "assembly.fasta").write_text('>synthetic\nACGT\n')
        for scenario in ["samples", "summary", "database", "other", "mob_new", "mob_existing", "mob_incomplete"]:
            case = tmp / scenario
            case.mkdir()
            db = case / "database"
            db.mkdir()
            if scenario == "mob_existing":
                (db / "clusters.txt").write_text('mock clusters\n')
                (db / "status.txt").write_text('complete\n')
            command = [nextflow, "-C", str(config), "run", str(tmp / "check.nf"),
                 "--scenario", scenario, "--outdir", str(case / "results"), "-ansi-log", "false",
                 "--analysis_type", "assembly" if scenario == "other" else "isolate-analysis",
                 "--mob_db", str(db), "--assembly", str(tmp / "assembly.fasta")]
            env = os.environ | {"NXF_OFFLINE": "true", "SOMATEM_MOB_INCOMPLETE": "1" if scenario == "mob_incomplete" else "0"}
            run = subprocess.run(
                command, cwd=case, env=env,
                text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            )
            assert (run.returncode == 0) == (scenario in ["samples", "mob_new", "mob_existing"]), run.stdout
            with (case / "trace.txt").open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            if scenario == "samples":
                assert set((case / "results/completed.txt").read_text().splitlines()) == {
                    "healthy:1", "transient:2"
                }, run.stdout
                assert sorted((r["tag"], r["exit"], r["attempt"]) for r in rows) == sorted([
                    ("healthy", "0", "1"), ("permanent", "1", "1"),
                    ("transient", "137", "1"), ("transient", "0", "2"),
                    ("exhausted", "137", "1"), ("exhausted", "137", "2"),
                ]), rows
                log = (case / ".nextflow.log").read_text()
                assert "executor 'local' > cpus=2; memory=1 GB;" in log, log
            elif scenario in ["mob_new", "mob_existing"]:
                assert len(rows) == 5 and all(r["exit"] == "0" for r in rows), rows
                assert len(list((case / "results").glob('isolate_analysis/*/plasmids/*.contig_report.txt'))) == 4
                if scenario == "mob_existing":
                    assert not (db / "init_calls.txt").exists()
                else:
                    assert (db / "init_calls.txt").read_text().splitlines() == ['init']
                    (db / "status.txt").unlink()
                    resumed = subprocess.run(command + ['-resume'], cwd=case, env=env,
                                             text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
                    assert resumed.returncode == 0, resumed.stdout
                    assert (db / "init_calls.txt").read_text().splitlines() == ['init', 'init']
            else:
                assert len(rows) == 1 and rows[0]["exit"] == "1", rows
                if scenario == "other":
                    assert "executor 'local' > cpus=2; memory=1 GB;" not in (case / ".nextflow.log").read_text()
            print(f"PASS {scenario}")


if __name__ == "__main__":
    main()
