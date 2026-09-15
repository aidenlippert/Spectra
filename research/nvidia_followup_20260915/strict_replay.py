"""Public replay caller with the even-N assumption of the spin proof enforced."""
import argparse
import json
from pathlib import Path


def require_supported_sector(data):
    m,n=data.get('modes'),data.get('particles')
    if type(m) is not int or type(n) is not int or m%2 or n%2 or not 2<=n<=m-2:
        raise ValueError('This complete spin-screen proof requires paired modes and even N, with a nonempty nonsinglet sector')


def run(case,proposal,output,rotated=None,compiled_grams=False,compiled_rationals=False):
    require_supported_sector(json.loads((case/'fixture.json').read_text()))
    from research.nvidia_followup_20260915.replay import run as replay
    return replay(case,proposal,output,rotated,compiled_grams,compiled_rationals)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('case',type=Path)
    p.add_argument('proposal')
    p.add_argument('output',type=Path)
    p.add_argument('--rotated',type=Path)
    p.add_argument('--compiled-grams',action='store_true')
    p.add_argument('--compiled-rationals',action='store_true')
    a=p.parse_args()
    run(a.case.resolve(),a.proposal,a.output.resolve(),a.rotated.resolve() if a.rotated else None,a.compiled_grams,a.compiled_rationals)
