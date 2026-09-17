"""Campaign-local ledger, reusing the existing bounded process runner."""
import argparse
from pathlib import Path
from research.correlated_pair_20260913 import budget

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--name',required=True);p.add_argument('--seconds',type=int,required=True)
    p.add_argument('command',nargs=argparse.REMAINDER);a=p.parse_args()
    budget.OUT=Path(__file__).resolve().parents[2]/'results/intervention_reduction_20260916/runs'
    command=a.command[1:] if a.command and a.command[0]=='--' else a.command
    raise SystemExit(budget.run(a.name,a.seconds,command))
