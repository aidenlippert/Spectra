"""Bind the measured library experiment to its preserved parent campaign."""
from datetime import datetime, timezone
from fractions import Fraction
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT/'results/nvidia_followup_20260915'
PARENT = ROOT/'results/transfer_solver_20260915'


def read(path):
    return json.loads(path.read_text())


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        while block := stream.read(1 << 20):
            digest.update(block)
    return digest.hexdigest()


def write(path, value):
    with path.open('x') as stream:
        json.dump(value, stream, indent=2)
        stream.write('\n')


def run():
    parent = read(PARENT/'seal.json')
    if not parent['all_hashes_match'] or sha(PARENT/'manifest.json') != parent['manifest_sha256']:
        raise ValueError('Parent campaign seal is not valid')
    stopped = read(OUT/'termination_confirmed.json')
    if stopped['status'] not in ('absent', 'terminated') or stopped['id'] != read(OUT/'instance.json')['id']:
        raise ValueError('Owned cloud instance termination is not confirmed')
    source = read(OUT/'source_upload_manifest.json')
    if sha(OUT/'source.tar.gz') != source['sha256']:
        raise ValueError('Uploaded source archive changed')
    downloaded = OUT/'download_final'
    if sha(OUT/'final_download.tar.gz') != read(OUT/'final_download_receipt.json')['archive_sha256']:
        raise ValueError('Downloaded result archive changed')
    download_files = read(downloaded/'final_download_manifest.json')
    for name, expected in download_files.items():
        if sha(downloaded/name) != expected:
            raise ValueError(('Downloaded result changed', name))
    protocol = read(downloaded/'results/nvidia_followup_20260915/accelerated_protocol.json')
    for name, expected in protocol['source_sha256'].items():
        path = ROOT/name
        if path.name == 'replay.py':
            path = OUT/'code_snapshots/replay_mixed_arithmetic_failed.py'
        if sha(path) != expected:
            raise ValueError(('Executed comparison source changed', name))
    retry = read(OUT/'replay_retry_upload.json')
    for name, expected in retry['files'].items():
        if sha(ROOT/'research/nvidia_followup_20260915'/name) != expected:
            raise ValueError(('Executed replay retry source changed', name))
    summary = read(OUT/'summary.json')
    for candidate in summary['all_H10_certificates']:
        interval = read(Path(candidate['path'])/'interval.json')
        width = Fraction(interval['upper_Ha'])-Fraction(interval['lower_Ha'])
        if not 0 <= width <= Fraction(16, 10000):
            raise ValueError('H10 exact width misses its target')
        if any(Fraction(interval[k]) != Fraction(candidate[k]) for k in ('upper_Ha', 'lower_Ha')):
            raise ValueError('Reported H10 endpoints changed')
    write(OUT/'final_validation.json', {
        'parent_manifest_sha256': parent['manifest_sha256'],
        'previous_files_verified_by_parent': parent['previous_files_verified'],
        'downloaded_files_verified': len(download_files),
        'executed_comparison_sources_verified': len(protocol['source_sha256']),
        'executed_retry_sources_verified': len(retry['files']),
        'H10_locally_verified_certificates': len(summary['all_H10_certificates']),
        'owned_cloud_instance_terminated': True,
        'scope': 'Final provenance and reported-endpoint checks; scientific replay is charged in the retained process receipts.'})
    files = {}
    for folder in (ROOT/'research/nvidia_followup_20260915', OUT):
        for path in sorted(folder.rglob('*')):
            if path.is_file() and '__pycache__' not in path.parts and path.name not in ('manifest.json', 'seal.json'):
                files[str(path.relative_to(ROOT))] = {'bytes': path.stat().st_size, 'sha256': sha(path)}
    write(OUT/'manifest.json', {
        'kind': 'measured_nvidia_library_followup_v1',
        'created_UTC': datetime.now(timezone.utc).isoformat(),
        'parent_manifest': str((PARENT/'manifest.json').relative_to(ROOT)),
        'parent_sha256': parent['manifest_sha256'],
        'files': files, 'file_count': len(files)})
    if any(sha(ROOT/name) != row['sha256'] for name, row in files.items()):
        raise ValueError('Experiment files changed while sealing')
    write(OUT/'seal.json', {'manifest_sha256': sha(OUT/'manifest.json'),
        'new_files_verified': len(files), 'all_hashes_match': True,
        'owned_instance_terminated': True})
    print(json.dumps(read(OUT/'seal.json')))


if __name__ == '__main__':
    run()
