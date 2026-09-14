"""Fresh stdlib-only replay of all stored training and transfer certificates."""
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import argparse
import hashlib
import json
import sys
import time

from research.side_routes_20260913.finite_range import replay


def run(directory, out, workers):
    start = time.monotonic()
    files = sorted((directory/'searches').glob('*.json')) + sorted((directory/'witnesses').glob('*.json'))
    payloads = []
    training = transfer = 0
    for path in files:
        record = json.loads(path.read_text())
        if 'candidates' in record:
            payloads.extend(record['candidates']); training += len(record['candidates'])
        else:
            payloads.append(record); transfer += 1
    if not payloads:
        raise ValueError('No exact claims to replay')
    with ProcessPoolExecutor(max_workers=workers) as pool:
        results = list(pool.map(replay,payloads))
    receipt = {'accepted':len(results),'training_candidates':training,'transfer_witnesses':transfer,
               'python':sys.version,'no_site':bool(sys.flags.no_site),'workers':workers,
               'wall_seconds':time.monotonic()-start,
               'sum_verification_seconds':sum(r['wall_seconds'] for r in results),
               'verifier_sha256':hashlib.sha256(Path(__file__).with_name('finite_range.py').read_bytes()).hexdigest(),
               'input_sha256':{str(p.relative_to(directory)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
               'scope':'Fresh execution of the same exact verifier. Independent CAR checks are separate focused tests.'}
    out.write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:receipt[k] for k in ('accepted','no_site','wall_seconds')}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--directory',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--workers',type=int,default=4)
    args = parser.parse_args()
    if not 1 <= args.workers <= 8: raise ValueError('Worker budget must be 1 through 8')
    run(args.directory,args.out,args.workers)
