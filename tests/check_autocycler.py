#!/usr/bin/env python3
"""Run with `python3 tests/check_autocycler.py [case ...]`; requires Nextflow, Java and Python.

Executes the isolate workflow and real module shell scripts
with fake scientific executables. Checks routing/contracts, not assembly quality.
"""

import json
import csv
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
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

if args in [['--version'], ['-V']] or (tool in ['metaMDBG', 'bwa'] and not args):
    print('Version: 1.3.1' if tool == 'metaMDBG' else tool + ' 1.0')
    sys.exit(0)

read_type = os.environ['SOMATEM_TEST_READ_TYPE']
failure = os.environ.get('SOMATEM_TEST_FAILURE', '')
if tool == 'checkm2':
    inputs = [Path(p) for p in args[args.index('--input') + 1:args.index('--output-directory')]]
    output = Path(value('--output-directory'))
    candidates = output.name == 'candidate_quality'
    assert len(inputs) == (int(os.environ['SOMATEM_TEST_CANDIDATES']) if candidates else 1), inputs
    assert len({p.name.split('.')[0] for p in inputs}) == 1, inputs
    assert value('--database_path') == 'checkm2.dmnd', args
    output.mkdir()
    rows = ['Name\tCompleteness\tContamination\n']
    for p in sorted(inputs, reverse=True):
        quality = {'flye': (99, 5), 'myloasm': (98, 1), 'metamdbg': (99.5, 3),
                   'canu': (93, 0), 'hifiasm': (97, 1)}[p.name.split('.')[1]] if candidates else (98, 1)
        if failure == 'invalid_quality':
            quality = ('nan', 0)
        rows.append(f'{p.stem}\t{quality[0]}\t{quality[1]}\n')
    if failure == 'missing_quality':
        rows.pop()
    (output / 'quality_report.tsv').write_text(''.join(rows))
    sys.exit(0)
if tool == 'kraken2':
    Path(value('--report')).write_text('100.00\t1\t1\tS\t1\tmock\n')
    sys.exit(0)
if tool == 'bakta':
    shutil.copyfile(args[-1], value('--prefix') + '.fna')
    sys.exit(0)
if tool == 'mob_recon':
    output = Path(value('--outdir'))
    output.mkdir()
    shutil.copyfile(value('--infile'), output / 'contig_report.txt')
    sys.exit(0)
if tool == 'btyper3':
    output = Path(value('-o')) / 'btyper3_final_results'
    output.mkdir(parents=True)
    shutil.copyfile(value('-i'), output / (Path(value('-i')).stem + '.txt'))
    sys.exit(0)
if tool == 'bwa':
    assert args[0] in ['index', 'mem'], args
    print('mock alignment')
    sys.exit(0)
if tool == 'polypolish':
    if args[0] == 'filter':
        shutil.copyfile(value('--in1'), value('--out1'))
        shutil.copyfile(value('--in2'), value('--out2'))
    else:
        assert args[0] == 'polish', args
        print(Path(args[1]).read_text(), end='')
    sys.exit(0)
if tool == 'pypolca':
    output = Path(value('-o'))
    output.mkdir()
    shutil.copyfile(value('-a'), output / 'pypolca_corrected.fasta')
    sys.exit(0)
if tool == 'autocycler' and args[0] in ['compress', 'cluster', 'trim', 'resolve', 'combine']:
    sample = next(Path('assemblies').glob('*.fasta')).name.split('.')[0]
    affected = sample == 'alpha' or failure == 'both'
    if affected and failure == args[0]:
        sys.exit('Simulated Autocycler command failure')
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
    fasta.write_text(f'>{sample}:{assembler}:{subset}:{read_type}\n' + 'A' * 25 + '\n')
    graph.write_text('S\t1\tACGT\n')
elif args[0] == 'subsample':
    output = Path(value('--out_dir'))
    output.mkdir()
    for i in range(1, int(value('--count')) + 1):
        shutil.copyfile(value('--reads'), output / f'sample_{i:02d}.fastq')
elif args[0] == 'compress':
    assert value('--max_contigs') == os.environ['SOMATEM_TEST_MAX_CONTIGS'], args
    inputs = sorted(Path(value('-i')).glob('*.fasta'))
    assert len(inputs) == int(os.environ['SOMATEM_TEST_CANDIDATES']), inputs
    headers = [p.read_text().splitlines()[0] for p in inputs]
    assert len(set(headers)) == len(inputs), headers
    assert len({h.split(':')[0] for h in headers}) == 1, headers
    output = Path(value('-a'))
    output.mkdir()
    (output / 'consensus_assembly.fasta').write_text(''.join(p.read_text() for p in inputs))
    if affected and failure:
        (output / 'consensus_assembly.fasta').write_text('>incomplete\nACGT\n')
    if affected and failure == 'boundary':
        (output / 'consensus_assembly.fasta').write_text('>boundary\n' + 'A' * 80 + '\n')
    if affected and failure == 'empty_consensus':
        (output / 'consensus_assembly.fasta').write_text('')
    (output / 'consensus_assembly.gfa').write_text('S\t1\tACGT\n')
elif args[0] == 'cluster':
    assert value('--max_contigs') == os.environ['SOMATEM_TEST_MAX_CONTIGS'], args
    if affected and failure == 'no_clusters':
        sys.exit(0)
    cluster = Path(value('-a')) / 'clustering/qc_pass/cluster_001'
    cluster.mkdir(parents=True)
    (cluster / '5_final.gfa').write_text('S\t1\tACGT\n')
else:
    assert args[0] in ['trim', 'resolve', 'combine'], args
'''


def main():
    nextflow = shutil.which("nextflow")
    assert nextflow, "Nextflow must be on PATH"
    defaults = "\n".join(re.findall(
        r"^    (?:autocycler_\w+|flye_mode|hybrid_assembly|run_\w+|polypolish_\w+|pypolca_\w+|min_shortread_depth_for_pypolca)\s*=.*$",
        (ROOT / "nextflow.config").read_text(), re.MULTILINE
    ))
    # Fail if the documented default list and the runtime configuration diverge.
    schema = json.loads((ROOT / "nextflow_schema.json").read_text())
    assert schema["properties"]["autocycler_assemblers"]["default"] == ",".join(ASSEMBLERS)
    max_contigs = schema["properties"]["autocycler_max_contigs"]["default"]
    assert f'autocycler_max_contigs     = {max_contigs}' in defaults

    with tempfile.TemporaryDirectory(prefix="somatem-autocycler-check-") as tmp:
        tmp = Path(tmp)
        binaries = tmp / "bin"
        binaries.mkdir()
        (binaries / 'python').symlink_to(sys.executable)
        (binaries / 'somatem_report.py').symlink_to(ROOT / 'bin/somatem_report.py')
        for tool in ["autocycler", "flye", "myloasm", "metaMDBG", "canu", "hifiasm", "checkm2", "bwa", "polypolish", "pypolca", "kraken2", "pigz", "bakta", "mob_recon", "btyper3"]:
            executable = binaries / tool
            executable.write_text(MOCK)
            executable.chmod(0o755)
        for sample in ["alpha", "beta"]:
            for suffix in ['', '_R1', '_R2']:
                (tmp / f"{sample}{suffix}.fastq").write_text(sample + "\n")
        (tmp / 'checkm2.dmnd').touch()
        (tmp / 'db').mkdir()
        (tmp / 'db/clusters.txt').write_text('mock database\n')
        (tmp / 'samples.csv').write_text('sample,long_reads\nalpha,alpha.fastq\nbeta,beta.fastq\n')
        workflow = tmp / "check.nf"
        workflow.write_text(f"include {{ ISOLATE_ANALYSIS }} from '{ROOT}/subworkflows/local/isolate_analysis'\n" + """
workflow {
    ISOLATE_ANALYSIS(channel.of(
        tuple([id: 'alpha'], file('alpha.fastq'), file('alpha_R1.fastq'), file('alpha_R2.fastq'), '100', ''),
        tuple([id: 'beta'], file('beta.fastq'), file('beta_R1.fastq'), file('beta_R2.fastq'), '100', '')
    ), channel.value(file('db')), channel.value(file('db')), channel.of(file('checkm2.dmnd')))
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
params.analysis_type = 'isolate-analysis'
params.output_dir = params.outdir
params.input = 'samples.csv'
params.mob_suite_db = "${projectDir}/db"
params.kraken2_save_output_fastqs = false
params.kraken2_save_reads_assignment = false
params.btyper3_args = ''
params.data_type = 'isolate'
params.sample_environment = ''
params.sequencing_technology = 'nanopore'
""" + "\n".join(
            f"params.autocycler_{name}_args = '--test-arg {name}'"
            for name in ASSEMBLERS if name != "flye"
        ))

        cases = [
            ("default", ASSEMBLERS, "ont_r10", [], None),
            ("max_contigs", ["flye"], "ont_r10",
             ["--autocycler_assemblers", "flye", "--autocycler_max_contigs", "400"], None),
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
        cases += [(name, ASSEMBLERS, 'ont_r10', [], None) for name in [
            'incomplete', 'empty_consensus', 'no_clusters', 'compress', 'cluster', 'trim', 'resolve', 'combine', 'both', 'boundary'
        ]]
        cases += [(name, ASSEMBLERS, 'ont_r10', options, None) for name, options in [
            ('hybrid', ['--hybrid_assembly', 'true']),
            ('pypolca_only', ['--hybrid_assembly', 'true', '--run_polypolish', 'false']),
            ('polypolish_only', ['--hybrid_assembly', 'true', '--run_pypolca', 'false']),
            ('no_polishing', ['--hybrid_assembly', 'true', '--run_polypolish', 'false', '--run_pypolca', 'false']),
            ('typing', ['--run_btyper3', 'true']),
        ]]
        cases += [(name, ASSEMBLERS, 'ont_r10', [], 'Error executing process')
                  for name in ['missing_quality', 'invalid_quality']]
        cases += [
            (f"bad_max_contigs_{i}", [], "ont_r10", ["--autocycler_max_contigs", value],
             "--autocycler_max_contigs must be a positive integer")
            for i, value in enumerate(["0", "-1", "1.5", "typo", "4294967296"])
        ]
        cases += [
            (f'bad_fraction_{i}', [], 'ont_r10', ['--autocycler_min_genome_fraction', value],
             '--autocycler_min_genome_fraction must be greater than 0 and at most 1')
            for i, value in enumerate(['0', '-1', '1.1', 'typo', 'NaN'])
        ]
        for name, assemblers, read_type, options, error in cases:
            if sys.argv[1:] and name not in sys.argv[1:]:
                continue
            output = tmp / name
            failure = ('incomplete' if name in ['hybrid', 'pypolca_only', 'polypolish_only', 'no_polishing', 'typing']
                       else name if name in ['incomplete', 'empty_consensus', 'no_clusters', 'compress', 'cluster', 'trim', 'resolve', 'combine', 'both', 'boundary', 'missing_quality', 'invalid_quality'] else '')
            env = dict(os.environ, PATH=f"{binaries}{os.pathsep}{os.environ['PATH']}",
                       NXF_OFFLINE="true", NXF_ANSI_LOG="false",
                       SOMATEM_TEST_READ_TYPE=read_type,
                       SOMATEM_TEST_FLYE_MODE="--nano-corr" if name == "legacy_corrected" else {
                           "ont_r10": "--nano-hq", "ont_r9": "--nano-raw",
                           "pacbio_hifi": "--pacbio-hifi", "pacbio_clr": "--pacbio-raw"
                       }[read_type],
                       SOMATEM_TEST_CANDIDATES=str(4 * len(assemblers)),
                       SOMATEM_TEST_MAX_CONTIGS="400" if name == "max_contigs" else str(max_contigs),
                       SOMATEM_TEST_FAILURE=failure,
                       SOMATEM_TEST_EMPTY="myloasm" if name == "empty_assembly" else "")
            result = subprocess.run(
                [nextflow, "-C", str(config), "run", str(workflow),
                 "-work-dir", str(tmp / f"work-{name}"), "--outdir", str(output), "--output_dir", str(output), *options],
                cwd=tmp, env=env, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                timeout=180
            )
            if error:
                assert result.returncode != 0 and error in result.stdout, result.stdout
            else:
                assert result.returncode == 0, result.stdout
                assert (output / 'summary_reports/isolate_analysis_summary_report.html').exists()
                for sample in ["alpha", "beta"]:
                    sample_dir = output / f'isolate_analysis/{sample}'
                    assert (sample_dir / f'autocycler/{sample}.autocycler.log').exists()
                    fasta = sample_dir / f'assembly/{sample}.fasta'
                    headers = {line for line in fasta.read_text().splitlines() if line.startswith(">")}
                    fallback = failure and failure != 'boundary' and (sample == 'alpha' or failure == 'both')
                    expected = {
                        f">{sample}:{assembler}:subsample_{i:02d}:{read_type}"
                        for assembler in assemblers for i in range(1, 5)
                    }
                    if fallback:
                        expected = {f'>{sample}:myloasm:subsample_01:{read_type}'}
                        selection = sample_dir / f'quality/assembly_selection/{sample}.assembly_selection.tsv'
                        with selection.open() as handle:
                            rows = list(csv.DictReader(handle, delimiter='\t'))
                        assert (selection.parent / f'{sample}.candidate_quality.tsv').exists()
                        assert len(rows) == 4 * len(assemblers), rows
                        assert [r['Assembly'] for r in rows if r['Selected'] == 'True'] == [f'{sample}.myloasm.subsample_01.assembly.fasta'], rows
                        assert not (sample_dir / f'assembly/{sample}.autocycler.fasta').exists()
                        assert not (sample_dir / f'assembly/{sample}.autocycler.gfa').exists()
                        assert (sample_dir / f'autocycler/{sample}.autocycler.warning.txt').exists()
                        assert f'Sample {sample}: Autocycler failed' in result.stdout, result.stdout
                        assert f'Sample {sample}: continuing with' in result.stdout, result.stdout
                    else:
                        assert not (sample_dir / 'quality/assembly_selection').exists()
                        if failure == 'boundary' and sample == 'alpha':
                            expected = {'>boundary'}
                    assert headers == expected, (name, sample, headers, expected)
                    assert (sample_dir / f'annotation/bakta/{sample}.fna').read_text() == fasta.read_text()
                    assert (sample_dir / f'plasmids/{sample}.contig_report.txt').read_text() == fasta.read_text()
                    assert (sample_dir / f'quality/checkm2/{sample}_checkm2_report.tsv').exists()
                    if name == 'typing':
                        assert (sample_dir / f'typing/btyper3/{sample}.txt').read_text() == fasta.read_text()
            print(f"PASS {name}", flush=True)


if __name__ == "__main__":
    main()
