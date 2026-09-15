"""Record bounded local attempts without modifying any sealed earlier campaign."""
import argparse
import fcntl
from pathlib import Path
import sys
from research.correlated_pair_20260913 import budget as runner
from research.transfer_solver_20260915.budget import dump

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT/'results/interacting_scaling_20260915'
STD = '/opt/homebrew/Caskroom/miniconda/base/bin/python'
NUM = str(ROOT/'.venv-correlated/bin/python')
CHEM = str(ROOT/'.venv-molecule/bin/python')


def run(name, seconds, command):
    if (OUT/'manifest.json').exists():
        raise RuntimeError('Preserve the sealed campaign')
    OUT.mkdir(exist_ok=True)
    with (OUT/'local_compute.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        runner.OUT = OUT/'runs'
        return runner.run(name, seconds, command)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--name', required=True)
    p.add_argument('--seconds', type=int, required=True)
    p.add_argument('command', nargs=argparse.REMAINDER)
    a = p.parse_args()
    command = a.command[1:] if a.command and a.command[0] == '--' else a.command
    sys.exit(run(a.name, a.seconds, command))
