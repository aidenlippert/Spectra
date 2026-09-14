"""Fresh output destinations for the existing bounded process runner."""
import argparse
from pathlib import Path
import sys
from research.correlated_pair_20260913 import budget as runner
runner.OUT=Path(__file__).resolve().parents[2]/'results/collective_completion_20260914/runs'
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--name',required=True);p.add_argument('--seconds',type=int,required=True);p.add_argument('command',nargs=argparse.REMAINDER);a=p.parse_args()
    sys.exit(runner.run(a.name,a.seconds,a.command[1:] if a.command and a.command[0]=='--' else a.command))
