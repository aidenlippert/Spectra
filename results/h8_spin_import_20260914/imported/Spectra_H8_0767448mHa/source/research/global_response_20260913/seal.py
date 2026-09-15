"""Preservation and broad local evidence inventory, not a compact proof format."""
from pathlib import Path
import argparse
import hashlib
import json
import time

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT/'results/global_response_20260913'


def digest(path):
    raw = path.read_bytes()
    return {'bytes': len(raw), 'sha256': hashlib.sha256(raw).hexdigest()}


def verify(files):
    for name, expected in files.items():
        path = ROOT/name
        if not path.is_file() or digest(path) != expected:
            raise ValueError('Evidence file changed or missing: '+name)


def run(create=False):
    start = time.monotonic(); path = OUT/'manifest.json'
    if create:
        baseline = json.loads((OUT/'preservation_before.json').read_text())['files']
        verify(baseline)
        files = dict(baseline)
        for folder in (ROOT/'research/global_response_20260913', OUT):
            for f in sorted(folder.rglob('*')):
                if f.is_file() and f != path and '__pycache__' not in f.parts and f.suffix != '.pyc' and f.name != '.DS_Store':
                    files[str(f.relative_to(ROOT))] = digest(f)
        receipt = {'prior_manifest_sha256': digest(ROOT/'results/composable_response_20260913/manifest.json')['sha256'],
                   'preserved_baseline_count': len(baseline), 'new_campaign_artifact_count': len(files)-len(baseline),
                   'files': files,
                   'scope': 'Broad preservation/evidence inventory, including prior unrelated research. Its size is not a compact mathematical certificate size.'}
        path.write_text(json.dumps(receipt, indent=2, sort_keys=True)+'\n')
    receipt = json.loads(path.read_text()); verify(receipt['files'])
    print(json.dumps({'verified_files': len(receipt['files']), 'preserved_baseline_count': receipt['preserved_baseline_count'],
                      'new_artifact_count': receipt['new_campaign_artifact_count'], 'manifest_sha256': digest(path)['sha256'],
                      'seconds': time.monotonic()-start}))


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('--create', action='store_true'); run(p.parse_args().create)
