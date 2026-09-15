"""Experimental FLINT rational backend; existing stdlib endpoints authoritative."""
import argparse
from fractions import Fraction
import importlib
import json
from pathlib import Path
import sys
import time


def run(case, output, stdlib=False, grams=False):
    if not stdlib or grams:
        import flint
    names=['experiments.marginal_symbolic','research.collective_completion_20260914.spin_replay',
        'research.collective_completion_20260914.spin_screen','research.certificate_scaling.spin_twirl',
        'research.correlated_pair_20260913.mps_exact']
    modules=[importlib.import_module(n) for n in names]
    # All scientific module text remains unchanged; this proposal experiment
    # compares the arithmetic implementation in a separate process.
    changed=[]
    for name,module in ([] if stdlib or grams else list(sys.modules.items())):
        if name.startswith(('research.','experiments.')):
            for alias in ('F','Fraction'):
                if getattr(module,alias,None) is Fraction:
                    setattr(module,alias,flint.fmpq)
                    changed.append(name+'.'+alias)
    if grams:
        from research.nvidia_followup_20260915.flint_squares import expand_squares
        modules[0].expand_squares=expand_squares
        changed.append('experiments.marginal_symbolic.expand_squares')
    read=lambda p:json.loads(p.read_text())
    expected=read(case/'interval.json')
    start=time.monotonic()
    lower=modules[2].check(read(case/'fixture.json'),read(case/'certificate.json'),read(case/'nonsinglet.json'))
    lower_seconds=time.monotonic()-start
    upper=modules[4].check(read(case/'fixture.json'),read(case/'mps/state.json'))
    if Fraction(lower['lower'])!=Fraction(expected['lower_Ha']) or Fraction(upper['upper_Ha'])!=Fraction(expected['upper_Ha']):
        raise ValueError('FLINT endpoints differ from the authoritative stdlib receipt')
    record={'library':'FLINT_integer_Gram' if grams else ('stdlib' if stdlib else 'python-flint'),'version':sys.version if stdlib and not grams else flint.__version__,'substituted_aliases':changed,
        'lower_seconds':lower_seconds,'complete_replay_seconds':time.monotonic()-start,
        'exact_endpoints_equal_to_stdlib':True,'experimental_backend_only':True}
    with output.open('x') as stream:
        json.dump(record,stream,indent=2)
    print(json.dumps(record),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('case',type=Path)
    p.add_argument('output',type=Path)
    p.add_argument('--stdlib',action='store_true')
    p.add_argument('--grams',action='store_true')
    a=p.parse_args()
    run(a.case.resolve(),a.output.resolve(),a.stdlib,a.grams)
