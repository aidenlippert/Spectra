"""Seal the complete local campaign and verify the preserved manifest chain."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
from research.interacting_scaling_20260915.budget import ROOT, OUT, dump

PARENT = 'results/nvidia_followup_20260915/manifest.json'
EXPECTED_PARENT = '34c01e6e5b112f90cb18b38fcc544b7dd99465cf20b251adb57d74177ded7599'


def sha(path):
    value = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024*1024), b''): value.update(chunk)
    return value.hexdigest()


def checked_path(relative):
    path = (ROOT/relative).resolve()
    if not path.is_relative_to(ROOT): raise ValueError('Manifest path escapes the workspace')
    return path


def verify_chain(relative, expected, visited=None):
    visited = set() if visited is None else visited
    if relative in visited: raise ValueError('Cyclic manifest ancestry')
    visited.add(relative)
    path = checked_path(relative)
    if sha(path) != expected: raise ValueError('Changed manifest: '+relative)
    manifest = json.loads(path.read_text())
    count = 0
    for name, record in manifest['files'].items():
        file = checked_path(name)
        digest = record['sha256'] if isinstance(record, dict) else record
        if sha(file) != digest: raise ValueError('Changed preserved file: '+name)
        if isinstance(record, dict) and 'bytes' in record and file.stat().st_size != record['bytes']:
            raise ValueError('Changed byte count: '+name)
        count += 1
    if manifest.get('parent_manifest'):
        count += verify_chain(manifest['parent_manifest'], manifest['parent_sha256'], visited)
    return count


def verify_preserved():
    count = verify_chain(PARENT, EXPECTED_PARENT)
    inventory = 'results/h8_spin_import_20260914/preservation_before.json'
    parent = json.loads(checked_path('results/h8_spin_import_20260914/manifest.json').read_text())
    if inventory not in parent['files']:
        raise ValueError('Historical preservation inventory is not bound by the checked parent chain')
    protected = json.loads(checked_path(inventory).read_text())['files']
    for name, record in protected.items():
        if sha(checked_path(name)) != record['sha256']:
            raise ValueError('Changed earlier preserved input: '+name)
    return count+len(protected)


def seal():
    if (OUT/'manifest.json').exists(): raise ValueError('This campaign is already sealed')
    running = [p.name for p in (OUT/'runs').glob('*.json') if json.loads(p.read_text())['status'] == 'starting']
    if running: raise ValueError('Cannot seal unfinished runs: '+str(running))
    preserved = verify_preserved()
    files = {}
    for root in (ROOT/'research/interacting_scaling_20260915', OUT):
        for path in sorted(root.rglob('*')):
            if path.is_file() and '__pycache__' not in path.parts and path.name not in ('manifest.json', 'seal.json', 'local_compute.lock'):
                files[str(path.relative_to(ROOT))] = {'sha256': sha(path), 'bytes': path.stat().st_size}
    dump(OUT/'manifest.json', {'kind': 'direct_interacting_certificate_campaign_v1',
        'created_UTC': datetime.now(timezone.utc).isoformat(), 'parent_manifest': PARENT,
        'parent_sha256': EXPECTED_PARENT, 'files': files, 'file_count': len(files)})
    digest = sha(OUT/'manifest.json')
    checked = verify_chain(str((OUT/'manifest.json').relative_to(ROOT)), digest)
    dump(OUT/'seal.json', {'manifest_sha256': digest, 'new_files_verified': len(files),
        'preserved_chain_files_verified': preserved, 'total_file_checks': checked,
        'file_integrity_is_separate_from_scientific_acceptance': True})
    print(json.dumps(json.loads((OUT/'seal.json').read_text()), indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('action', choices=('seal', 'check', 'preserved'))
    a = p.parse_args()
    if a.action == 'seal': seal()
    elif a.action == 'preserved': print(verify_preserved())
    else:
        receipt = json.loads((OUT/'seal.json').read_text())
        print(verify_chain(str((OUT/'manifest.json').relative_to(ROOT)), receipt['manifest_sha256']))
