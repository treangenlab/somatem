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
import time

tool = Path(sys.argv[0]).name
args = sys.argv[1:]
def value(flag):
    return args[args.index(flag) + 1]

if args in [['--version'], ['-V']] or (tool in ['metaMDBG', 'bwa'] and not args):
    print('Version: 1.3.1' if tool == 'metaMDBG' else tool + ' 1.0')
    sys.exit(0)

read_type = os.environ['SOMATEM_TEST_READ_TYPE']
failure = os.environ.get('SOMATEM_TEST_FAILURE', '')
length_case = os.environ.get('SOMATEM_TEST_LENGTH_CASE', '')
if tool == 'checkm2':
    inputs = [Path(p) for p in args[args.index('--input') + 1:args.index('--output-directory')]]
    output = Path(value('--output-directory'))
    assert len(inputs) == 1, inputs
    assert value('--database_path') == 'checkm2.dmnd', args
    output.mkdir()
    (output / 'quality_report.tsv').write_text(
        f'Name\tCompleteness\tContamination\n{inputs[0].stem}\t98\t1\n'
    )
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
        expected_size = 'invalid' if failure == 'invalid_genome_size' and Path(value('--reads')).read_text().strip() == 'alpha' else '100'
        assert value('--genome_size') == expected_size, args
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
    marker = os.environ.get('SOMATEM_TEST_CONSENSUS_MARKER')
    if marker and sample == 'beta' and assembler == 'flye' and subset == 'subsample_01':
        # Leave the other executor slots free; alpha must progress before beta finishes.
        deadline = time.monotonic() + 30
        while not Path(marker).exists():
            assert time.monotonic() < deadline, 'Alpha consensus waited for unfinished beta assembly'
            time.sleep(0.1)
    length = 100
    if sample == 'alpha' and length_case:
        if length_case == 'all_short_candidates' or assembler != 'flye':
            length = 79
        elif subset == 'subsample_02':
            length = 80
    header = f'>{sample}:{assembler}:{subset}:{read_type}'
    sequence = header + '\n' + 'A' * length + '\n'
    if length_case == 'short_candidates' and sample == 'alpha' and assembler == 'flye' and subset == 'subsample_01':
        # Both contigs are short individually; their whole-assembly length passes.
        sequence = header + '\n' + 'A' * 40 + '\n' + header + ':part2\n' + 'A' * 60 + '\n'
    fasta.write_text(sequence)
    graph.write_text('S\t1\tACGT\n')
elif args[0] == 'subsample':
    output = Path(value('--out_dir'))
    output.mkdir()
    for i in range(1, int(value('--count')) + 1):
        shutil.copyfile(value('--reads'), output / f'sample_{i:02d}.fastq')
    if failure == 'invalid_genome_size' and Path(value('--reads')).read_text().strip() == 'alpha':
        Path('genome_size.txt').write_text('invalid\n')
elif args[0] == 'compress':
    assert value('--threads') == '1', args
    assert value('--max_contigs') == os.environ['SOMATEM_TEST_MAX_CONTIGS'], args
    inputs = sorted(Path(value('-i')).glob('*.fasta'))
    expected_count = 4 if length_case == 'short_candidates' and sample == 'alpha' else int(os.environ['SOMATEM_TEST_CANDIDATES'])
    assert len(inputs) == expected_count, inputs
    if length_case == 'all_short_candidates' and sample == 'alpha':
        sys.exit('Short candidates reached compress')
    headers = [p.read_text().splitlines()[0] for p in inputs]
    assert len(set(headers)) == len(inputs), headers
    assert len({h.split(':')[0] for h in headers}) == 1, headers
    marker = os.environ.get('SOMATEM_TEST_CONSENSUS_MARKER')
    if marker and sample == 'alpha':
        Path(marker).touch()
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
    minimum = os.environ['SOMATEM_TEST_MIN_ASSEMBLIES']
    if minimum:
        assert value('--min_assemblies') == minimum, args
    else:
        assert '--min_assemblies' not in args, args
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
    ), channel.value(file('db')), channel.value(file('db')), channel.value(file('checkm2.dmnd')))
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
            ("short_candidates", ASSEMBLERS, "ont_r10", [], None),
            ("two_subsamples", ASSEMBLERS, "ont_r10", ["--autocycler_subsamples", "2"], None),
            ("single_assembly", ["flye"], "ont_r10",
             ["--autocycler_assemblers", "flye", "--autocycler_subsamples", "1"], None),
            ("min_assemblies", ASSEMBLERS, "ont_r10", ["--autocycler_min_assemblies", "3"], None),
            ("min_assemblies_strict", ASSEMBLERS, "ont_r10", ["--autocycler_min_assemblies", "5"], None),
            ("streaming", ["flye", "canu"], "ont_r10",
             ["--autocycler_assemblers", "flye,canu", "--autocycler_subsamples", "2"], None),
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
        consensus_failures = {
            'all_short_candidates': 'No candidate assemblies meet the minimum genome fraction',
            'invalid_genome_size': 'Invalid',
            'incomplete': 'Incomplete Autocycler consensus',
            'empty_consensus': 'Incomplete Autocycler consensus',
            'no_clusters': 'No Autocycler QC-pass clusters',
            'both': 'Incomplete Autocycler consensus',
            **{name: 'Simulated Autocycler command failure' for name in ['compress', 'cluster', 'trim', 'resolve', 'combine']},
        }
        cases += [(name, ASSEMBLERS, 'ont_r10', [], error) for name, error in consensus_failures.items()]
        cases += [('boundary', ASSEMBLERS, 'ont_r10', [], None)]
        cases += [(name, ASSEMBLERS, 'ont_r10', options, None) for name, options in [
            ('hybrid', ['--hybrid_assembly', 'true']),
            ('pypolca_only', ['--hybrid_assembly', 'true', '--run_polypolish', 'false']),
            ('polypolish_only', ['--hybrid_assembly', 'true', '--run_pypolca', 'false']),
            ('no_polishing', ['--hybrid_assembly', 'true', '--run_polypolish', 'false', '--run_pypolca', 'false']),
            ('typing', ['--run_btyper3', 'true']),
        ]]
        cases += [
            (f"bad_max_contigs_{i}", [], "ont_r10", ["--autocycler_max_contigs", value],
             "--autocycler_max_contigs must be a positive integer")
            for i, value in enumerate(["0", "-1", "1.5", "typo", "4294967296"])
        ]
        cases += [
            (f"bad_min_assemblies_{i}", [], "ont_r10", ["--autocycler_min_assemblies", value],
             "--autocycler_min_assemblies must be a positive integer")
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
            subsamples = int(options[options.index('--autocycler_subsamples') + 1]) if '--autocycler_subsamples' in options else 4
            min_assemblies = {'min_assemblies': '3', 'min_assemblies_strict': '5'}.get(name, '')
            failure = name if name in consensus_failures or name == 'boundary' else ''
            env = dict(os.environ, PATH=f"{binaries}{os.pathsep}{os.environ['PATH']}",
                       NXF_OFFLINE="true", NXF_ANSI_LOG="false",
                       SOMATEM_TEST_READ_TYPE=read_type,
                       SOMATEM_TEST_FLYE_MODE="--nano-corr" if name == "legacy_corrected" else {
                           "ont_r10": "--nano-hq", "ont_r9": "--nano-raw",
                           "pacbio_hifi": "--pacbio-hifi", "pacbio_clr": "--pacbio-raw"
                       }[read_type],
                       SOMATEM_TEST_CANDIDATES=str(subsamples * len(assemblers)),
                       SOMATEM_TEST_CONSENSUS_MARKER=str(tmp / 'alpha-consensus-started') if name == 'streaming' else '',
                       SOMATEM_TEST_MAX_CONTIGS="400" if name == "max_contigs" else str(max_contigs),
                       SOMATEM_TEST_MIN_ASSEMBLIES=str(min_assemblies),
                       SOMATEM_TEST_LENGTH_CASE=name if name in ['short_candidates', 'all_short_candidates'] else '',
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
                if name in consensus_failures:
                    assert "Error executing process > 'ISOLATE_ANALYSIS:AUTOCYCLER_CONSENSUS" in result.stdout, result.stdout
                    for sample in ['alpha', 'beta'] if name == 'both' else ['alpha']:
                        sample_dir = output / f'isolate_analysis/{sample}'
                        assert not list((sample_dir / 'assembly').glob('*')), sample_dir
                        assert not (sample_dir / 'annotation').exists(), sample_dir
                        assert not (sample_dir / 'quality/checkm2').exists(), sample_dir
                    if name == 'all_short_candidates':
                        audits = list((tmp / f'work-{name}').glob('*/*/alpha.autocycler.lengths.tsv'))
                        assert len(audits) == 1, audits
                        with audits[0].open() as handle:
                            lengths = list(csv.DictReader(handle, delimiter='\t'))
                        assert len(lengths) == subsamples * len(assemblers), lengths
                        assert all(row['Stage'] == 'candidate' and row['Passes_length_filter'] == 'False' for row in lengths), lengths
            else:
                assert result.returncode == 0, result.stdout
                assert (output / 'summary_reports/isolate_analysis_summary_report.html').exists()
                for sample in ["alpha", "beta"]:
                    sample_dir = output / f'isolate_analysis/{sample}'
                    assert (sample_dir / f'autocycler/{sample}.autocycler.log').exists()
                    fasta = sample_dir / f'assembly/{sample}.fasta'
                    headers = {line for line in fasta.read_text().splitlines() if line.startswith(">")}
                    expected = {
                        f">{sample}:{assembler}:subsample_{i:02d}:{read_type}"
                        for assembler in assemblers for i in range(1, subsamples + 1)
                    }
                    if name == 'short_candidates' and sample == 'alpha':
                        expected = {f'>{sample}:flye:subsample_{i:02d}:{read_type}' for i in range(1, 5)}
                        expected.add(f'>{sample}:flye:subsample_01:{read_type}:part2')
                    with (sample_dir / f'autocycler/{sample}.autocycler.lengths.tsv').open() as handle:
                        lengths = list(csv.DictReader(handle, delimiter='\t'))
                    candidates = [row for row in lengths if row['Stage'] == 'candidate']
                    assert len(candidates) == subsamples * len(assemblers), candidates
                    assert all(int(row['Expected_bp']) == 100 for row in lengths), lengths
                    assert all(float(row['Genome_fraction']) == int(row['Length_bp']) / 100 for row in lengths), lengths
                    for row in candidates:
                        original = sample_dir / 'autocycler/assemblies' / Path(row['Assembly']).name
                        assert sum(len(line.strip()) for line in original.read_text().splitlines() if not line.startswith('>')) == int(row['Length_bp']), original
                    if name == 'short_candidates' and sample == 'alpha':
                        assert sum(row['Passes_length_filter'] == 'True' for row in candidates) == 4, candidates
                        assert sorted(int(row['Length_bp']) for row in candidates) == [79] * 16 + [80, 100, 100, 100], candidates
                    else:
                        assert all(row['Passes_length_filter'] == 'True' for row in candidates), candidates
                    consensus = [row for row in lengths if row['Stage'] == 'consensus']
                    assert len(consensus) == 1 and consensus[0]['Passes_length_filter'] == 'True', consensus
                    assert int(consensus[0]['Length_bp']) == sum(len(line.strip()) for line in fasta.read_text().splitlines() if not line.startswith('>')), consensus
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
