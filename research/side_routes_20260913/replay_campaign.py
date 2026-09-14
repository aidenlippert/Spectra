"""Replay every saved exact chain result using the stdlib-only verifier."""
import argparse
import hashlib
import json
import sys
import time
from pathlib import Path
from research.side_routes_20260913.positive_chain import replay


def run(directory, out):
    start = time.monotonic()
    screen_path = directory / 'screen.json'
    files = sorted(directory.glob('*_m*.json'))
    inputs = [(str(screen_path), item) for item in json.loads(screen_path.read_text())]
    inputs.extend((str(path), json.loads(path.read_text())) for path in files)
    for _, payload in inputs:
        replay(payload)
    receipt = {'accepted': len(inputs), 'screen_count': len(inputs) - len(files),
               'held_out_size_count': len(files), 'python': sys.version, 'no_site': bool(sys.flags.no_site),
               'verifier_sha256': hashlib.sha256(Path(__file__).with_name('positive_chain.py').read_bytes()).hexdigest(),
               'input_sha256': {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in [screen_path] + files},
               'wall_seconds': time.monotonic() - start,
               'scope': 'Fresh replay of the same exact algorithm. Independent CAR enumeration is in the focused tests, not this replay.'}
    out.write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({k: receipt[k] for k in ('accepted', 'no_site', 'wall_seconds')}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--directory', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    run(args.directory, args.out)
