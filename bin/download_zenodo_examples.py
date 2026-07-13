#!/usr/bin/env python3
"""Download and arrange Somatem example datasets from a Zenodo record."""

import argparse
import hashlib
import shutil
import sys
import tempfile
import urllib.request
import zipfile
from pathlib import Path


DATASETS = (
    ('16S.zip', '16S', '586e03b8f84bc8f44b21d85a8a6f01ee'),
    ('mag_assembly.zip', 'assembly', '93b1252c8b787841c63e6ca07f6e14e3'),
    ('isolate_assembly.zip', 'isolate_assembly', '055dd2612e91b0cba8daa9d03942d475'),
    ('metagenome_small.zip', 'metagenome_small', 'cb6b8af82c927f6b691b651d144ce234'),
    ('other_tools_files.zip', 'other_tools_files', '865f1368b1cd8f289dae2ca314326d49'),
    ('seqscreen.zip', 'other_tools_files/more_files/seqscreen', '014db703b23fbcd04e1e1d8c9b28ed9b'),
)

SAMPLESHEETS = {
    '16S_pilot_samples.csv': (
        'sample,fastq_1',
        ('mock1', '16S/mockm95_sub10k.fastq.gz'),
        ('mock2', '16S/mockm91_sub10k.fastq.gz'),
    ),
    'meta_tax_samples.csv': (
        'sample,fastq_1',
        ('mock9', 'metagenome_small/mock9_sub10k.fastq.gz'),
        ('mock20', 'metagenome_small/mock20_sub10k.fastq.gz'),
    ),
    'mag_samples.csv': (
        'sample,fastq_1',
        ('zymo', 'assembly/mock20_hiq100k.fastq.gz'),
    ),
    'isolate_samples.csv': (
        'sample,long_reads,short_reads_1,short_reads_2,expected_genome_size,species',
        (
            'FFI_BCgr',
            'isolate_assembly/FFI_BCgr36_ONT_q15_min1000.fastq.gz',
            'isolate_assembly/FFI_BCgr36_Ill_1.fastq.gz',
            'isolate_assembly/FFI_BCgr36_Ill_2.fastq.gz',
            '5.4m',
            'Bacillus_cereus',
        ),
    ),
    'timeseries_samples.csv': (
        'sample,fastq_1',
        ('t0', 'other_tools_files/rhea/t0.fasta'),
        ('t1', 'other_tools_files/rhea/t1.fasta'),
    ),
    'seqscreen_samples.csv': (
        'sample,fasta',
        ('example', 'other_tools_files/more_files/seqscreen/10_seqs.fasta'),
    ),
}


def md5sum(path):
    digest = hashlib.md5()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def download(url, destination):
    request = urllib.request.Request(
        url, headers={'User-Agent': 'somatem-example-data-downloader'}
    )
    with urllib.request.urlopen(request) as response, destination.open('wb') as handle:
        shutil.copyfileobj(response, handle)


def copy_contents(source, destination):
    destination.mkdir(parents=True, exist_ok=True)
    for item in source.iterdir():
        target = destination / item.name
        if item.is_dir():
            shutil.copytree(item, target, dirs_exist_ok=True)
        else:
            shutil.copy2(item, target)


def unpack(archive, destination, work_dir):
    extracted = work_dir / archive.stem
    with zipfile.ZipFile(archive) as handle:
        handle.extractall(extracted)

    entries = list(extracted.iterdir())
    source = entries[0] if len(entries) == 1 and entries[0].is_dir() else extracted
    copy_contents(source, destination)


def write_samplesheets(save_dir):
    data_dir = (save_dir / 'data').resolve()
    samplesheet_dir = save_dir / 'samplesheets'
    samplesheet_dir.mkdir(parents=True, exist_ok=True)

    for filename, (header, *records) in SAMPLESHEETS.items():
        lines = [header]
        for record in records:
            fields = [
                str(data_dir / value) if '/' in value else value
                for value in record
            ]
            lines.append(','.join(fields))
        (samplesheet_dir / filename).write_text('\n'.join(lines) + '\n')


def file_url(record_id, filename):
    return f'https://zenodo.org/api/records/{record_id}/files/{filename}/content'


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--record-id', required=True)
    parser.add_argument('--save-dir', required=True, type=Path)
    parser.add_argument('--manifest', required=True, type=Path)
    parser.add_argument('--samplesheets-only', action='store_true')
    args = parser.parse_args()

    args.save_dir = args.save_dir.resolve()
    if args.samplesheets_only:
        write_samplesheets(args.save_dir)
        return

    data_dir = args.save_dir / 'data'
    completed = []

    with tempfile.TemporaryDirectory(prefix='somatem-zenodo-') as temporary:
        temporary_dir = Path(temporary)
        for filename, target_name, expected_md5 in DATASETS:
            target_dir = data_dir / target_name
            if target_dir.is_dir() and any(target_dir.iterdir()):
                print(f'Reusing existing example data: {target_dir}', flush=True)
                completed.append(str(target_dir))
                continue

            archive = temporary_dir / filename
            print(f'Downloading {filename} from Zenodo...', flush=True)
            download(file_url(args.record_id, filename), archive)

            print(f'Verifying {filename}...', flush=True)
            actual_md5 = md5sum(archive)
            if actual_md5 != expected_md5:
                raise RuntimeError(f'Checksum failed for {filename}: {actual_md5}')

            print(f'Extracting {filename} to {target_dir}...', flush=True)
            unpack(archive, target_dir, temporary_dir)
            completed.append(str(target_dir))

    write_samplesheets(args.save_dir)
    args.manifest.write_text('\n'.join(completed) + '\n')


if __name__ == '__main__':
    try:
        main()
    except Exception as error:
        print(f'ERROR: {error}', file=sys.stderr)
        sys.exit(1)
