"""One bounded process at a time; preserve all successful and failed attempts."""
import argparse
import json
from pathlib import Path
import sys
from research.correlated_pair_20260913 import budget as runner

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT/'results/transfer_solver_20260915'


def dump(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x') as stream:
        json.dump(value, stream, indent=2)
        stream.write('\n')


if __name__ == '__main__':
    if (OUT/'manifest.json').exists():
        raise RuntimeError('Sealed completed pass')
    p = argparse.ArgumentParser()
    p.add_argument('--name', required=True)
    p.add_argument('--seconds', type=int, required=True)
    p.add_argument('command', nargs=argparse.REMAINDER)
    a = p.parse_args()
    runner.OUT = OUT/'runs'
    cmd = a.command[1:] if a.command and a.command[0] == '--' else a.command
    sys.exit(runner.run(a.name, a.seconds, cmd))
