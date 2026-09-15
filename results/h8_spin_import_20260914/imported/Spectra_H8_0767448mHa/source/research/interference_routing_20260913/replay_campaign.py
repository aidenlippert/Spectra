"""Fresh-process replay of every saved interval, using standard library only."""
from pathlib import Path
import argparse
import hashlib
import json
import time

from research.interference_routing_20260913.routing import replay as molecular_replay
from research.interference_routing_20260913.local_templates import replay as local_replay
from research.interference_routing_20260913.sector_counterexample import run as check_counterexample


def run(directory):
    started = time.monotonic()
    molecular = json.loads((directory/'molecular/summary.json').read_text())
    local = json.loads((directory/'local/summary.json').read_text())
    records = []
    for item in molecular['cases']:
        if item['status'] != 'accepted':
            continue
        path = directory/'molecular'/(item['name']+'_proof.json')
        actual = molecular_replay(item['fixture'], json.loads(path.read_text()))
        for key in ('lower', 'upper', 'width'):
            if actual[key] != item['receipt'][key]:
                raise ValueError('Molecular receipt mismatch')
        records.append({'proof': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
                        'width': actual['width'], 'seconds': actual['replay_seconds']})
    for item in local['cases']:
        path = directory/'local'/item['file']
        actual = local_replay(json.loads(path.read_text()))
        for key in ('lower', 'upper', 'width'):
            if actual[key] != item['receipt'][key]:
                raise ValueError('Local receipt mismatch')
        records.append({'proof': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
                        'width': actual['width'], 'seconds': actual['replay_seconds']})
    counterexample = json.loads((directory/'sector_counterexample.json').read_text())
    if counterexample != json.loads(json.dumps(check_counterexample())):
        raise ValueError('Sector counterexample receipt mismatch')
    return {'status': 'all_replayed', 'count': len(records), 'records': records,
            'sector_counterexample_replayed': True,
            'wall_seconds': time.monotonic()-started,
            'scope': 'Fresh process, same accepting algorithms; independent CAR checks are in test_routing.py.'}


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--directory', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    args = p.parse_args()
    result = run(args.directory)
    args.out.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'records'}))
