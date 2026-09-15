"""Bounded, separately recorded checks of the supplied H8 certificate."""
import argparse
from pathlib import Path
import sys
from research.correlated_pair_20260913 import budget as runner

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT/'results/h8_spin_import_20260914'

if __name__ == '__main__':
    if (OUT/'manifest.json').exists():
        raise RuntimeError('Sealed validation pass')
    p = argparse.ArgumentParser()
    p.add_argument('--name', required=True)
    p.add_argument('--seconds', required=True, type=int)
    p.add_argument('command', nargs=argparse.REMAINDER)
    args = p.parse_args()
    command = args.command[1:] if args.command and args.command[0] == '--' else args.command
    runner.OUT = OUT/'runs'
    sys.exit(runner.run(args.name, args.seconds, command))
