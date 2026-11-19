"""Run with Nextflow, Java and DIAMOND on PATH. No SeqScreen database needed."""

import os
from pathlib import Path
import re
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[1]


def main():
    with tempfile.TemporaryDirectory(prefix="somatem-seqscreen-") as directory:
        tmp = Path(directory)
        db = tmp / "db" / "diamond"
        db.mkdir(parents=True)
        (tmp / "reference.faa").write_text(">reference\n" + "MKLYEAVRNPQWSTDHGICF" * 4 + "\n")
        codons = "ATG AAA CTG TAC GAA GCT GTT CGT AAC CCT CAA TGG TCT ACT GAT CAT GGT ATT TGT TTT"
        (tmp / "query.fa").write_text(">query\n" + codons.replace(" ", "") * 4 + "\n")
        subprocess.run(["diamond", "makedb", "--in", "reference.faa", "-d", str(db / "uniref.mini")],
                       cwd=tmp, check=True)
        (tmp / "check.config").write_text(f"""
params.analysis_type = 'seqscreen'
params.threads = 2
params.max_memory_gb = 1
params.local_memory_reserve_gb = 4
params.seqscreen_evalue = 30
includeConfig '{ROOT}/conf/base.config'
""")
        (tmp / "check.nf").write_text(f"""
include {{ SEQSCREEN_DIAMOND_BLASTX }} from '{ROOT}/modules/local/seqscreen/diamond_blastx/main'
workflow {{
    samples = channel.of('a', 'b').map {{ id -> tuple([id: id], file('query.fa')) }}
    SEQSCREEN_DIAMOND_BLASTX(samples, file('db'))
}}
""")
        subprocess.run(["nextflow", "-C", "check.config", "run", "check.nf"], cwd=tmp, check=True,
                       env=os.environ | {"NXF_OFFLINE": "true"})
        log = (tmp / ".nextflow.log").read_text()
        assert "executor 'local' > cpus=2; memory=1 GB;" in log
        assert "peakRunning=1;" in log  # Each sample requests the entire memory budget.
        for sample in ("a", "b"):
            (tabular,) = (tmp / "work").glob(f"*/*/{sample}.ur100.btab")
            assert tabular.read_text().startswith("query\treference\t")
            assert "<Hit>" in tabular.with_suffix(".xml").read_text()
            script = (tabular.parent / ".command.sh").read_text()
            assert re.findall(r"--block-size (\S+)", script) == ["0.125", "0.125"]
            assert "--index-chunks 1" not in script
        print("PASS SeqScreen local budget and both DIAMOND outputs for two samples")


if __name__ == "__main__":
    main()
