#!/usr/bin/env python3
"""Run with `python3 tests/check_autocycler.py`; requires Nextflow, Java and Python.

Executes the isolate workflow's assembly section and real module shell scripts
with fake scientific executables. Checks routing/contracts, not assembly quality.
"""

import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[1]
ASSEMBLERS = ["flye", "myloasm", "metamdbg", "canu", "hifiasm"]
MOCK = r'''#!/usr/bin/env python3
import os
from pathlib import Path
import shutil
import sys

tool = Path(sys.argv[0]).name
args = sys.argv[1:]
def value(flag):
    return args[args.index(flag) + 1]

if args == ['--version'] or (tool == 'metaMDBG' and not args):
    print('Version: 1.3.1' if tool == 'metaMDBG' else tool + ' 1.0')
    sys.exit(0)

read_type = os.environ['SOMATEM_TEST_READ_TYPE']
if tool == 'flye' or args[0] == 'helper':
    assembler = tool if tool == 'flye' else args[1]
    if tool == 'flye':
        mode = os.environ['SOMATEM_TEST_FLYE_MODE']
        reads = Path(value(mode))
        output = Path(value('--out-dir'))
        output.mkdir()
        fasta = output / 'assembly.fasta'
        graph = output / 'assembly_graph.gfa'
        (output / 'flye.log').touch()
        (output / 'assembly_info.txt').touch()
    else:
        assert value('--read_type') == read_type, args
        assert value('--genome_size') == '100', args
        assert value('--test-arg') == assembler, args
        if assembler == 'canu':
            assert 'maxMemory=1' in args, args
        reads = Path(value('--reads'))
        prefix = value('--out_prefix')
        fasta = Path(prefix + '.fasta')
        graph = Path(prefix + '.gfa')
        if os.environ.get('SOMATEM_TEST_EMPTY') == assembler:
            fasta.touch()
            sys.exit(0)
    assert value('--threads') == '1', args
    sample = reads.read_text().strip()
    subset = reads.stem.replace('sample_', 'subsample_')
    fasta.write_text(f'>{sample}:{assembler}:{subset}:{read_type}\nACGT\n')
    graph.write_text('S\t1\tACGT\n')
elif args[0] == 'subsample':
    output = Path(value('--out_dir'))
    output.mkdir()
    for i in range(1, int(value('--count')) + 1):
        shutil.copyfile(value('--reads'), output / f'sample_{i:02d}.fastq')
elif args[0] == 'compress':
    inputs = sorted(Path(value('-i')).glob('*.fasta'))
    assert len(inputs) == int(os.environ['SOMATEM_TEST_CANDIDATES']), inputs
    headers = [p.read_text().splitlines()[0] for p in inputs]
    assert len(set(headers)) == len(inputs), headers
    assert len({h.split(':')[0] for h in headers}) == 1, headers
    output = Path(value('-a'))
    output.mkdir()
    (output / 'consensus_assembly.fasta').write_text(''.join(p.read_text() for p in inputs))
    (output / 'consensus_assembly.gfa').write_text('S\t1\tACGT\n')
elif args[0] == 'cluster':
    cluster = Path(value('-a')) / 'clustering/qc_pass/cluster_001'
    cluster.mkdir(parents=True)
    (cluster / '5_final.gfa').write_text('S\t1\tACGT\n')
else:
    assert args[0] in ['trim', 'resolve', 'combine'], args
'''


def main():
    nextflow = shutil.which("nextflow")
    assert nextflow, "Nextflow must be on PATH"
    source = (ROOT / "subworkflows/local/isolate_analysis.nf").read_text()
    # Exercise the actual assembly statements without unrelated database/annotation jobs.
    includes = "\n".join(
        line.replace("../../modules/", f"{ROOT}/modules/")
        for line in source.splitlines()
        if line.startswith("include {") and "AUTOCYCLER" in line
    )
    body = source.split("    main:\n", 1)[1].split("    def run_polypolish", 1)[0]
    kraken = body.index("    ch_kraken_reads =")
    body = body[:kraken] + body[body.index("    ch_autocycler_input ="):]
    defaults = "\n".join(re.findall(
        r"^    (?:autocycler_\w+|flye_mode)\s*=.*$",
        (ROOT / "nextflow.config").read_text(), re.MULTILINE
    ))
    # Fail if the documented default list and the runtime configuration diverge.
    schema = json.loads((ROOT / "nextflow_schema.json").read_text())
    assert schema["properties"]["autocycler_assemblers"]["default"] == ",".join(ASSEMBLERS)

    with tempfile.TemporaryDirectory(prefix="somatem-autocycler-check-") as tmp:
        tmp = Path(tmp)
        binaries = tmp / "bin"
        binaries.mkdir()
        for tool in ["autocycler", "flye", "myloasm", "metaMDBG", "canu", "hifiasm"]:
            executable = binaries / tool
            executable.write_text(MOCK)
            executable.chmod(0o755)
        for sample in ["alpha", "beta"]:
            (tmp / f"{sample}.fastq").write_text(sample + "\n")
        workflow = tmp / "check.nf"
        workflow.write_text(includes + "\nworkflow ASSEMBLY_CHECK {\n    take:\n    ch_reads\n    main:\n" + body + "\n}\n" + """
workflow {
    ASSEMBLY_CHECK(channel.of(
        tuple([id: 'alpha'], file('alpha.fastq'), '', '', '100', ''),
        tuple([id: 'beta'], file('beta.fastq'), '', '', '100', '')
    ))
}
""")
        config = tmp / "check.config"
        config.write_text("params {\n" + defaults + "\n    publish_dir_mode = 'copy'\n}\n" + """
process.cpus = 1
process.memory = '1 GB'
process.executor = 'local'
process.shell = ['bash', '-C', '-e', '-u', '-o', 'pipefail']
executor.queueSize = 4
conda.enabled = false
docker.enabled = false
singularity.enabled = false
""" + "\n".join(
            f"params.autocycler_{name}_args = '--test-arg {name}'"
            for name in ASSEMBLERS if name != "flye"
        ))

        cases = [
            ("default", ASSEMBLERS, "ont_r10", [], None),
            ("subset", ["canu", "hifiasm"], "pacbio_hifi",
             ["--autocycler_assemblers", " Canu, HIFIASM, canu ", "--autocycler_read_type", "pacbio_hifi"], None),
            ("legacy_hifi", ["flye"], "pacbio_hifi",
             ["--autocycler_assemblers", "flye", "--flye_mode", "pacbio-hifi"], None),
            ("explicit_hifi", ["flye", "myloasm", "metamdbg"], "pacbio_hifi",
             ["--autocycler_assemblers", "flye,myloasm,metamdbg", "--autocycler_read_type", "pacbio_hifi"], None),
            ("legacy_corrected", ["flye", "canu"], "ont_r9",
             ["--autocycler_assemblers", "flye,canu", "--flye_mode", "nano-corr"], None),
            ("unknown", [], "ont_r10", ["--autocycler_assemblers", "flye,typo"], "comma-separated selection"),
            ("empty", [], "ont_r10", ["--autocycler_assemblers", " , "], "comma-separated selection"),
            ("bad_read_type", [], "ont_r10", ["--autocycler_read_type", "typo"], "must be ont_r10"),
            ("low_accuracy", [], "ont_r9", ["--autocycler_read_type", "ont_r9"], "require high-accuracy reads"),
            ("empty_assembly", ["myloasm"], "ont_r10", ["--autocycler_assemblers", "myloasm"], "Error executing process"),
        ]
        for name, assemblers, read_type, options, error in cases:
            output = tmp / name
            env = dict(os.environ, PATH=f"{binaries}{os.pathsep}{os.environ['PATH']}",
                       NXF_OFFLINE="true", NXF_ANSI_LOG="false",
                       SOMATEM_TEST_READ_TYPE=read_type,
                       SOMATEM_TEST_FLYE_MODE="--nano-corr" if name == "legacy_corrected" else {
                           "ont_r10": "--nano-hq", "ont_r9": "--nano-raw",
                           "pacbio_hifi": "--pacbio-hifi", "pacbio_clr": "--pacbio-raw"
                       }[read_type],
                       SOMATEM_TEST_CANDIDATES=str(4 * len(assemblers)),
                       SOMATEM_TEST_EMPTY="myloasm" if name == "empty_assembly" else "")
            result = subprocess.run(
                [nextflow, "-C", str(config), "run", str(workflow),
                 "-work-dir", str(tmp / f"work-{name}"), "--outdir", str(output), *options],
                cwd=tmp, env=env, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                timeout=180
            )
            if error:
                assert result.returncode != 0 and error in result.stdout, result.stdout
            else:
                assert result.returncode == 0, result.stdout
                for sample in ["alpha", "beta"]:
                    fasta = output / f"isolate_analysis/{sample}/assembly/{sample}.autocycler.fasta"
                    headers = {line for line in fasta.read_text().splitlines() if line.startswith(">")}
                    assert headers == {
                        f">{sample}:{assembler}:subsample_{i:02d}:{read_type}"
                        for assembler in assemblers for i in range(1, 5)
                    }, headers
            print(f"PASS {name}", flush=True)


if __name__ == "__main__":
    main()
