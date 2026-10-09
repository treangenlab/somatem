#!/usr/bin/env python3
"""Download Somatem example data from Zenodo, quickly and reliably.

Version 2 of the somatem/bin/download_zenodo_examples.py.
This new version downloads distinct ZIP files concurrently and,
when the server supports byte ranges, downloads each large ZIP in parallel
chunks. Downloaded chunks persist across interrupted runs. Verified archives
are extracted concurrently, except for seqscreen.zip which must follow
other_tools_files.zip because their destination paths overlap.
"""

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import sys
import tempfile
import threading
import time
import urllib.error
import urllib.request
import zipfile

DATASETS = (
    ('16S.zip', '16S', '586e03b8f84bc8f44b21d85a8a6f01ee'),
    ('mag_assembly.zip', 'assembly', '93b1252c8b787841c63e6ca07f6e14e3'),
    ('isolate_assembly.zip', 'isolate_assembly', '055dd2612e91b0cba8daa9d03942d475'),
    ('metagenome_small.zip', 'metagenome_small', 'cb6b8af82c927f6b691b651d144ce234'),
    ('other_tools_files.zip', 'other_tools_files', '865f1368b1cd8f289dae2ca314326d49'),
    ('seqscreen.zip', 'other_tools_files/more_files/seqscreen', '014db703b23fbcd04e1e1d8c9b28ed9b'),
)
SAMPLESHEETS = {
    '16S_pilot_samples.csv': ('sample,fastq_1',
        ('mock1', '16S/mockm95_sub10k.fastq.gz'),
        ('mock2', '16S/mockm91_sub10k.fastq.gz')),
    'meta_tax_samples.csv': ('sample,fastq_1',
        ('mock9', 'metagenome_small/mock9_sub10k.fastq.gz'),
        ('mock20', 'metagenome_small/mock20_sub10k.fastq.gz')),
    'mag_samples.csv': ('sample,fastq_1',
        ('zymo', 'assembly/mock20_hiq100k.fastq.gz')),
    'isolate_samples.csv': ('sample,long_reads,short_reads_1,short_reads_2,expected_genome_size,species',
        ('FFI_BCgr', 'isolate_assembly/FFI_BCgr36_ONT_q15_min1000.fastq.gz',
         'isolate_assembly/FFI_BCgr36_Ill_1.fastq.gz',
         'isolate_assembly/FFI_BCgr36_Ill_2.fastq.gz', '5.4m', 'Bacillus_cereus')),
    'timeseries_samples.csv': ('sample,fastq_1',
        ('t0', 'other_tools_files/rhea/t0.fasta'),
        ('t1', 'other_tools_files/rhea/t1.fasta')),
    'seqscreen_samples.csv': ('sample,fasta',
        ('example', 'other_tools_files/more_files/seqscreen/10_seqs.fasta')),
}

CHUNK = 4 * 1024 * 1024
AGENT = 'somatem-example-data-downloader/3.0'
PRINT_LOCK = threading.Lock()
class RangeUnsupported(RuntimeError):
    """The server did not honor a requested HTTP byte range."""


RANGE_RE = re.compile(r'^bytes (\d+)-(\d+)/(\d+)$')


def log(message):
    with PRINT_LOCK:
        print(message, flush=True)


def md5sum(path):
    digest = hashlib.md5()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(CHUNK), b''):
            digest.update(block)
    return digest.hexdigest()


def file_url(record_id, filename):
    return f'https://zenodo.org/api/records/{record_id}/files/{filename}/content'


def open_url(url, timeout, headers=None):
    request = urllib.request.Request(
        url, headers={'User-Agent': AGENT, 'Accept-Encoding': 'identity', **(headers or {})}
    )
    return urllib.request.urlopen(request, timeout=timeout)


def write_samplesheets(save_dir):
    data_dir = (save_dir / 'data').resolve()
    samplesheet_dir = save_dir / 'samplesheets'
    samplesheet_dir.mkdir(parents=True, exist_ok=True)
    for filename, (header, *records) in SAMPLESHEETS.items():
        rows = [header]
        for record in records:
            rows.append(','.join(str(data_dir / value) if '/' in value else value for value in record))
        (samplesheet_dir / filename).write_text('\n'.join(rows) + '\n')


def retry_call(function, label, retries):
    for attempt in range(retries + 1):
        try:
            return function()
        except (OSError, RuntimeError, ValueError, urllib.error.URLError) as exc:
            if attempt >= retries:
                raise RuntimeError(f'{label} failed after {retries + 1} attempt(s): {exc}') from exc
            delay = min(20, 2 ** attempt)
            log(f'[retry] {label}: {exc}; sleeping {delay}s')
            time.sleep(delay)


def probe_ranges(url, timeout):
    """Return (total_bytes, etag) if Range is honored, otherwise None.

    No parallel segments are attempted without a strict 206 / Content-Range.
    """
    try:
        with open_url(url, timeout, {'Range': 'bytes=0-0'}) as response:
            if response.status != 206:
                return None
            match = RANGE_RE.fullmatch(response.headers.get('Content-Range', ''))
            if not match or int(match.group(1)) != 0 or int(match.group(2)) != 0:
                return None
            total = int(match.group(3))
            if total < 1:
                return None
            return total, response.headers.get('ETag', '')
    except (OSError, RuntimeError, urllib.error.URLError):
        return None


def download_sequential(url, archive, cache_dir, retries, timeout):
    partial = cache_dir / (archive.name + '.part')

    def attempt():
        offset = partial.stat().st_size if partial.exists() else 0
        headers = {'Range': f'bytes={offset}-'} if offset else {}
        try:
            response = open_url(url, timeout, headers)
        except urllib.error.HTTPError as error:
            if error.code == 416 and offset:
                # A completed .part is finalized by outer checksum verification.
                return partial
            raise
        with response:
            if offset and response.status == 206:
                match = RANGE_RE.fullmatch(response.headers.get('Content-Range', ''))
                if not match or int(match.group(1)) != offset:
                    raise RuntimeError('Server returned an incorrect resume range')
                mode = 'ab'
            elif response.status == 200:
                mode = 'wb'  # Server ignored Range: restart rather than append.
            else:
                raise RuntimeError(f'Unexpected HTTP {response.status}')
            expected = response.headers.get('Content-Length')
            received = 0
            with partial.open(mode) as output:
                while True:
                    block = response.read(CHUNK)
                    if not block:
                        break
                    output.write(block)
                    received += len(block)
            if expected is not None and received != int(expected):
                raise IOError(f'Incomplete response ({received} of {expected} bytes)')
        return partial

    log(f'[download] {archive.name}: single connection (resumable)')
    return retry_call(attempt, archive.name, retries)


def segment_download(url, part_path, start, end, total, etag, timeout, retries):
    expected = end - start + 1

    def attempt():
        already = part_path.stat().st_size if part_path.exists() else 0
        if already == expected:
            return part_path
        if already > expected:
            part_path.unlink()
            already = 0
        headers = {'Range': f'bytes={start + already}-{end}'}
        if etag and not etag.startswith('W/'):
            headers['If-Range'] = etag
        with open_url(url, timeout, headers) as response:
            if response.status != 206:
                raise RangeUnsupported('Remote server stopped supporting byte-range requests')
            match = RANGE_RE.fullmatch(response.headers.get('Content-Range', ''))
            if not match or tuple(map(int, match.groups())) != (start + already, end, total):
                raise RangeUnsupported('Remote server returned incorrect Content-Range')
            received = 0
            with part_path.open('ab' if already else 'wb') as output:
                while True:
                    block = response.read(CHUNK)
                    if not block:
                        break
                    output.write(block)
                    received += len(block)
            if part_path.stat().st_size != expected:
                raise IOError(f'Segment incomplete ({part_path.stat().st_size} of {expected} bytes)')
        return part_path

    # Range violations are not transient: immediately fall back to normal transfer.
    try:
        return retry_call(attempt, part_path.name, retries)
    except RuntimeError as exc:
        if isinstance(exc.__cause__, RangeUnsupported):
            raise RangeUnsupported(str(exc)) from exc
        raise


def download_segmented(url, archive, cache_dir, total, etag, parts, retries, timeout):
    pieces_dir = cache_dir / (archive.name + '.segments')
    pieces_dir.mkdir(parents=True, exist_ok=True)
    metadata_file = pieces_dir / 'metadata.json'
    new_meta = {'url': url, 'size': total, 'etag': etag, 'parts': parts}
    try:
        old_meta = json.loads(metadata_file.read_text())
    except (FileNotFoundError, ValueError):
        old_meta = None
    if old_meta != new_meta:
        # A changed record/version or layout makes previous pieces unsafe to reuse.
        shutil.rmtree(pieces_dir)
        pieces_dir.mkdir()
        metadata_file.write_text(json.dumps(new_meta))

    boundaries = [(total * i // parts, total * (i + 1) // parts - 1) for i in range(parts)]
    log(f'[download] {archive.name}: {parts} parallel range requests, {total / 1048576:.1f} MiB')
    with ThreadPoolExecutor(max_workers=parts) as executor:
        futures = [executor.submit(segment_download, url, pieces_dir / f'{i:03d}.part',
                                   start, end, total, etag, timeout, retries)
                   for i, (start, end) in enumerate(boundaries)]
        for future in as_completed(futures):
            future.result()

    assembled = cache_dir / (archive.name + '.assembling')
    try:
        with assembled.open('wb') as destination:
            for i in range(parts):
                with (pieces_dir / f'{i:03d}.part').open('rb') as source:
                    shutil.copyfileobj(source, destination, CHUNK)
        if assembled.stat().st_size != total:
            raise IOError(f'Merged archive is not {total} bytes')
        return assembled
    except Exception:
        assembled.unlink(missing_ok=True)
        raise


def download_one(record_id, dataset, cache_dir, parts, min_part_bytes, retries, timeout, force):
    filename, _, expected_md5 = dataset
    archive = cache_dir / filename
    url = file_url(record_id, filename)
    segments = cache_dir / (filename + '.segments')
    if force:
        archive.unlink(missing_ok=True)
        (cache_dir / (filename + '.part')).unlink(missing_ok=True)
        shutil.rmtree(segments, ignore_errors=True)

    if archive.is_file():
        if md5sum(archive) == expected_md5:
            log(f'[cached] {filename} verified')
            return archive
        log(f'[invalid] {filename}: removing bad cached archive')
        archive.unlink()

    old_partial = cache_dir / (filename + '.part')
    ranged = None if old_partial.exists() else probe_ranges(url, timeout)
    assembled = None
    if ranged and parts > 1:
        total, etag = ranged
        segment_count = min(parts, max(1, total // min_part_bytes))
        if segment_count > 1:
            try:
                assembled = download_segmented(url, archive, cache_dir, total, etag,
                                               segment_count, retries, timeout)
                if md5sum(assembled) != expected_md5:
                    raise IOError('MD5 checksum mismatch after assembling segments')
                os.replace(assembled, archive)
                shutil.rmtree(segments, ignore_errors=True)
            except (RangeUnsupported, IOError) as exc:
                # Only range incompatibility / corrupt merged data should trigger
                # single-connection fallback. Transient network failures retain
                # completed segments for the next run.
                if isinstance(exc, IOError) and not str(exc).startswith('MD5 checksum mismatch'):
                    raise
                log(f'[fallback] {filename}: {exc}; using one connection')
                if assembled:
                    assembled.unlink(missing_ok=True)
                shutil.rmtree(segments, ignore_errors=True)

    if not archive.exists():
        # Sequential retry uses existing partial file if one is present.
        partial = download_sequential(url, archive, cache_dir, retries, timeout)
        actual = md5sum(partial)
        if actual != expected_md5:
            partial.unlink(missing_ok=True)
            raise IOError(f'MD5 mismatch for {filename}: expected {expected_md5}, got {actual}')
        os.replace(partial, archive)
    log(f'[verified] {filename} ({archive.stat().st_size / 1048576:.1f} MiB)')
    return archive


def validate_zip(archive):
    for member in archive.infolist():
        path = PurePosixPath(member.filename)
        if (not member.filename or path.is_absolute() or '..' in path.parts or
                '\\' in member.filename or '\x00' in member.filename or
                (member.create_system == 3 and
                 (member.external_attr >> 16) & 0o170000 == 0o120000)):
            raise ValueError(f'Unsafe ZIP member: {member.filename!r}')


def unpack(archive, destination, data_dir):
    data_dir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.somatem-extract-', dir=data_dir) as tmp:
        extracted = Path(tmp) / 'contents'
        extracted.mkdir()
        with zipfile.ZipFile(archive) as zf:
            validate_zip(zf)
            zf.extractall(extracted)
        entries = list(extracted.iterdir())
        if not entries:
            raise ValueError(f'Empty archive: {archive}')
        source = entries[0] if len(entries) == 1 and entries[0].is_dir() else extracted
        destination.parent.mkdir(parents=True, exist_ok=True)
        if destination.exists():
            if destination.is_symlink() or not destination.is_dir():
                raise ValueError(f'Refusing to overwrite non-directory path: {destination}')
            shutil.copytree(source, destination, dirs_exist_ok=True)
        else:
            os.replace(source, destination)  # Fast rename on same filesystem.
    log(f'[extracted] {archive.name} -> {destination}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--record-id', required=True)
    parser.add_argument('--save-dir', required=True, type=Path)
    parser.add_argument('--manifest', required=True, type=Path)
    parser.add_argument('--samplesheets-only', action='store_true')
    parser.add_argument('--workers', type=int, default=3,
                        help='Concurrent archive downloads (default: 3)')
    parser.add_argument('--parts', type=int, default=3,
                        help='Concurrent HTTP byte-range requests per archive (default: 3)')
    parser.add_argument('--extract-workers', type=int, default=2,
                        help='Parallel ZIP extractions (default: 2)')
    parser.add_argument('--min-part-mb', type=int, default=32,
                        help='Minimum size of each HTTP range, MiB (default: 32)')
    parser.add_argument('--retries', type=int, default=4)
    parser.add_argument('--timeout', type=float, default=90)
    parser.add_argument('--cache-dir', type=Path,
                        help='Persistent download cache (default: SAVE_DIR/.somatem_download_cache)')
    parser.add_argument('--keep-archives', action='store_true')
    parser.add_argument('--force', action='store_true')
    args = parser.parse_args()
    if not (1 <= args.workers <= 16 and 1 <= args.parts <= 16 and
            1 <= args.extract_workers <= 8 and args.min_part_mb >= 1 and
            args.retries >= 0 and args.timeout > 0):
        parser.error('Invalid concurrency, size, retries, or timeout argument')

    save_dir = args.save_dir.expanduser().resolve()
    if args.samplesheets_only:
        write_samplesheets(save_dir)
        return
    data_dir = save_dir / 'data'
    cache_dir = (args.cache_dir or save_dir / '.somatem_download_cache').expanduser().resolve()
    data_dir.mkdir(parents=True, exist_ok=True)
    cache_dir.mkdir(parents=True, exist_ok=True)

    needed = []
    for dataset in DATASETS:
        _, target, _ = dataset
        destination = data_dir / target
        if destination.is_dir() and any(destination.iterdir()) and not args.force:
            log(f'[reuse] {destination} (not revalidated)')
        else:
            needed.append(dataset)

    # Extract completed files while others are still downloading.
    # seqscreen is intentionally deferred until after other_tools_files because
    # its output is nested under the other archive's destination.
    seqscreen_archive = None
    extraction_futures = {}
    errors = []
    with ThreadPoolExecutor(max_workers=args.extract_workers) as extract_pool:
        if needed:
            with ThreadPoolExecutor(max_workers=min(args.workers, len(needed))) as download_pool:
                futures = {
                    download_pool.submit(download_one, args.record_id, dataset, cache_dir,
                                         args.parts, args.min_part_mb * 1048576,
                                         args.retries, args.timeout, args.force): dataset
                    for dataset in needed
                }
                for future in as_completed(futures):
                    filename, target, _ = futures[future]
                    try:
                        archive = future.result()
                        if filename == 'seqscreen.zip':
                            seqscreen_archive = archive
                        else:
                            extraction_futures[extract_pool.submit(unpack, archive,
                                data_dir / target, data_dir)] = archive
                    except Exception as exc:
                        errors.append(f'{filename}: {exc}')

        for future in as_completed(extraction_futures):
            archive = extraction_futures[future]
            try:
                future.result()
                if not args.keep_archives:
                    archive.unlink(missing_ok=True)
            except Exception as exc:
                errors.append(f'Extraction {archive.name}: {exc}')

        if seqscreen_archive is not None and not errors:
            try:
                unpack(seqscreen_archive, data_dir / 'other_tools_files/more_files/seqscreen', data_dir)
                if not args.keep_archives:
                    seqscreen_archive.unlink(missing_ok=True)
            except Exception as exc:
                errors.append(f'Extraction seqscreen.zip: {exc}')

    if errors:
        raise RuntimeError('; '.join(errors))

    write_samplesheets(save_dir)
    manifest = args.manifest.expanduser()
    manifest.parent.mkdir(parents=True, exist_ok=True)
    manifest.write_text('\n'.join(str(data_dir / target) for _, target, _ in DATASETS) + '\n')
    log(f'Done: sample sheets in {save_dir / "samplesheets"}; manifest in {manifest}')


if __name__ == '__main__':
    try:
        main()
    except Exception as exc:
        print(f'ERROR: {exc}', file=sys.stderr)
        sys.exit(1)

