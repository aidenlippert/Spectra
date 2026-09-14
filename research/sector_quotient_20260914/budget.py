"""Durable bounded runs for this campaign; inherited output stays sealed."""
import argparse
from pathlib import Path
import sys
from research.correlated_pair_20260913 import budget as runner

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'results/sector_quotient_20260914'

if __name__ == '__main__':
    if (OUT / 'manifest.json').exists():
        raise RuntimeError('This campaign is sealed')
    p = argparse.ArgumentParser()
    p.add_argument('--name', required=True)
    p.add_argument('--seconds', required=True, type=int)
    p.add_argument('command', nargs=argparse.REMAINDER)
    a = p.parse_args()
    runner.OUT = OUT / 'runs'
    sys.exit(runner.run(a.name, a.seconds, a.command[1:] if a.command[0] == '--' else a.command))
