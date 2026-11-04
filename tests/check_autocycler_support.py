#!/usr/bin/env python3
"""Run with Python; requires Nextflow, Java and real Autocycler on PATH.

Checks candidate length filtering using deterministic random toy sequences only.
No reads, reference genomes, or other biological datasets are used.
"""

import csv
import os
from pathlib import Path
import random
import re
import shutil
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[1]


def main():
    nextflow = shutil.which("nextflow")
    assert nextflow and shutil.which("autocycler"), "Nextflow and Autocycler must be on PATH"
    defaults = "\n".join(re.findall(
        r"^    autocycler_(?:subsamples|min_assemblies|max_contigs|min_genome_fraction)\s*=.*$",
        (ROOT / "nextflow.config").read_text(), re.MULTILINE
    ))
    rng = random.Random(42)
    large = "".join(rng.choices("ACGT", k=50000))
    small = "".join(rng.choices("ACGT", k=2000))
    env = dict(os.environ, NXF_OFFLINE="true", NXF_ANSI_LOG="false")

    with tempfile.TemporaryDirectory(prefix="somatem-autocycler-support-") as tmp:
        for name, full_count, short_count, minimum, error in [
            ("three_full_seventeen_short", 3, 17, None, None),
            ("three_full_four_short", 3, 4, None, None),
            ("all_short", 0, 7, None, "No candidate assemblies meet the minimum genome fraction"),
            ("explicit_four", 3, 4, 4, "No Autocycler QC-pass clusters"),
        ]:
            case = Path(tmp) / name
            inputs = case / "inputs"
            inputs.mkdir(parents=True)
            for index in range(full_count + short_count):
                sequences = f">random_small\n{small}\n"
                if index < full_count:
                    sequences = f">random_large\n{large}\n" + sequences
                (inputs / f"candidate_{index:02}.fasta").write_text(sequences)
            (case / "genome_size.txt").write_text("52000\n")
            (case / "check.nf").write_text(
                f"include {{ AUTOCYCLER_CONSENSUS }} from '{ROOT}/modules/local/autocycler/consensus/main'\n"
                "workflow {\n"
                "    AUTOCYCLER_CONSENSUS(channel.of(tuple([id: 'synthetic'],\n"
                "        files('inputs/*.fasta'), file('genome_size.txt'))))\n"
                "}\n"
            )
            config = case / "check.config"
            config.write_text("params {\n" + defaults + "\n}\n" + """
params.outdir = 'results'
params.publish_dir_mode = 'copy'
process.cpus = 1
process.memory = '1 GB'
process.executor = 'local'
process.shell = ['bash', '-C', '-e', '-u', '-o', 'pipefail']
conda.enabled = false
docker.enabled = false
singularity.enabled = false
""")
            options = [] if minimum is None else ["--autocycler_min_assemblies", str(minimum)]
            run = subprocess.run(
                [nextflow, "-C", str(config), "run", "check.nf", *options],
                cwd=case, env=env, text=True, stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT, timeout=180,
            )
            published = case / "results/isolate_analysis/synthetic/assembly/synthetic.autocycler.fasta"
            if error is None:
                assert run.returncode == 0, run.stdout
            else:
                assert run.returncode != 0 and error in run.stdout, run.stdout
                assert not published.exists(), name
                assert not list((case / "work").glob("*/*/synthetic.autocycler.fasta")), name
            task_dirs = list((case / "work").glob("*/*/.command.sh"))
            assert len(task_dirs) == 1, run.stdout
            task = task_dirs[0].parent
            if not full_count:
                assert not (task / "autocycler_out").exists(), name
            staged = sorted((task / "assemblies").glob("*.fasta"))
            assert len(staged) == full_count, (name, staged)
            for path in staged:
                assert path.read_text() == (inputs / path.name).read_text(), path
            audit = task / "synthetic.autocycler.lengths.tsv"
            with audit.open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            candidates = [row for row in rows if row["Stage"] == "candidate"]
            assert len(candidates) == full_count + short_count, rows
            passing = {row["Assembly"] for row in candidates if row["Passes_length_filter"].lower() == "true"}
            assert passing == {f"candidate_{i:02}.fasta" for i in range(full_count)}, rows
            for row in candidates:
                assert int(row["Expected_bp"]) == 52000, row
                assert int(row["Length_bp"]) == (52000 if row["Assembly"] in passing else 2000), row
                assert (float(row["Genome_fraction"]) >= 0.8) == (row["Assembly"] in passing), row
            consensus = [row for row in rows if row["Stage"] == "consensus"]
            if error is None:
                lengths = sorted(len("".join(record.splitlines()[1:]))
                                 for record in published.read_text().split(">")[1:])
                assert lengths == [2000, 50000], (name, lengths)
                assert len(consensus) == 1 and int(consensus[0]["Length_bp"]) == 52000, rows
                assert consensus[0]["Passes_length_filter"].lower() == "true", rows
                published_audit = case / "results/isolate_analysis/synthetic/autocycler/synthetic.autocycler.lengths.tsv"
                assert published_audit.read_text() == audit.read_text()
            else:
                assert not consensus, rows
            print(f"PASS {name}: {full_count} complete candidates retained; "
                  f"{'52000-base consensus' if error is None else 'rejected'}", flush=True)


if __name__ == "__main__":
    main()
