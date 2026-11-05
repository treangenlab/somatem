#!/usr/bin/env python3
"""Run with `python3 tests/check_isolate_report.py`; no external tools or data needed."""

import csv
import html
import importlib.util
from pathlib import Path
import re
import subprocess
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("somatem_report", ROOT / "bin/somatem_report.py")
report = importlib.util.module_from_spec(spec)
spec.loader.exec_module(report)


def main():
    with tempfile.TemporaryDirectory(prefix="somatem-isolate-report-") as tmp:
        tmp = Path(tmp)
        inputs = tmp / "inputs"
        inputs.mkdir()
        samples = ["consensus", "fallback", "no_assembly", "candidate_only", "annotation_failed", "empty_final", "no_record"]
        samples += [f"missing_{index:03d}" for index in range(600)]
        samplesheet = tmp / "samples.csv"
        with samplesheet.open("w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(["sample", "long_reads"])
            writer.writerows([sample, "a" * 500 + ".fastq"] for sample in samples)
        assert samplesheet.stat().st_size > 200_000
        header, sample_rows = report.read_samplesheet(samplesheet)
        assert [row[0] for row in sample_rows] == samples

        for sample in samples[:7]:
            if sample != "no_record":
                with (inputs / f"{sample}.autocycler.selection.tsv").open("w", newline="") as handle:
                    writer = csv.writer(handle, delimiter="\t")
                    writer.writerow(["Sample", "Status", "Source", "Reason"])
                    status = sample if sample in {"fallback", "no_assembly"} else "consensus"
                    writer.writerow([sample, status, "candidate<script>.fasta", "consensus < minimum"])
            if sample in {"no_assembly", "candidate_only"}:
                (inputs / f"{sample}.autocycler.fasta").write_text(">candidate\nACGT\n")
                continue
            for suffix in [".fasta", ".kraken2.report.txt", "_checkm2_report.tsv", ".tsv", ".json", ".contig_report.txt"]:
                (inputs / f"{sample}{suffix}").write_text("{}\n" if suffix == ".json" else "output\n")
        (inputs / "annotation_failed.tsv").write_text("")
        (inputs / "empty_final.fasta").write_text("")

        # Outcome detection must not read assemblies, even for large cohorts.
        read_table_rows = report.read_table_rows
        def selection_rows(path, **kwargs):
            assert path.name.endswith(".autocycler.selection.tsv"), path
            return read_table_rows(path, **kwargs)
        report.read_table_rows = selection_rows
        section = report.isolate_status_section(list(inputs.iterdir()), header, sample_rows)
        rows = [
            [html.unescape(cell) for cell in re.findall(r"<td>(.*?)</td>", row)]
            for row in re.findall(r"<tr>(.*?)</tr>", section)
        ]
        outcomes = {row[0]: row for row in rows if row}
        assert len(outcomes) == len(samples), len(outcomes)
        assert outcomes["consensus"][1] == "Expected outputs present"
        assert outcomes["fallback"][1:3] == ["Expected outputs present", "fallback"]
        assert outcomes["no_assembly"][1] == "No usable assembly"
        assert outcomes["candidate_only"][1] == "Incomplete"
        assert "Final assembly" in outcomes["candidate_only"][-1]
        assert outcomes["annotation_failed"][1] == "Incomplete"
        assert outcomes["annotation_failed"][-1] == "Bakta table"
        assert outcomes["empty_final"][1] == "Incomplete"
        assert outcomes["no_record"][1] == "Incomplete"
        assert outcomes["missing_599"][1] == "Incomplete"
        assert "<script>" not in section and "candidate&lt;script&gt;.fasta" in section

        output = tmp / "report.html"
        subprocess.run([
            sys.executable, str(ROOT / "bin/somatem_report.py"),
            "--workflow", "isolate_analysis", "--workflow-label", "Isolate analysis",
            "--samplesheet", str(samplesheet), "--input-dir", str(inputs), "--output", str(output),
        ], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        document = output.read_text()
        assert '<section id="sample-outcomes">' in document
        assert "missing_599" in document
        assert "No usable assembly: 1" in document
    print("PASS isolate report: fallback provenance, missing/empty outputs, and all 607 input samples")


if __name__ == "__main__":
    main()
