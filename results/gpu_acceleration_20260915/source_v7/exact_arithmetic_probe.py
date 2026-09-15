"""Experimental GMP rational substitution; compare against stdlib replay.

No accepting-module files are edited. This isolated process substitutes only
their rational arithmetic alias. The independent stdlib result is authoritative.
"""
import argparse
from fractions import Fraction
import importlib
import json
from pathlib import Path
import time
import sys


def run(case, candidate, integers=False):
    import gmpy2
    names = ['experiments.marginal_symbolic',
             'research.collective_completion_20260914.spin_replay',
             'research.collective_completion_20260914.spin_screen',
             'research.certificate_scaling.spin_twirl',
             'research.correlated_pair_20260913.mps_exact']
    modules = [importlib.import_module(name) for name in names]
    substituted = []
    for name, module in list(sys.modules.items()):
        if name.startswith(('research.', 'experiments.')):
            for alias in ('F', 'Fraction'):
                if getattr(module, alias, None) is Fraction:
                    setattr(module, alias, gmpy2.mpq)
                    substituted.append(name+'.'+alias)
    if integers:
        original_step = modules[4].State.step
        def integer_step(self, environment, site, operator):
            # Validation still reads the original Python-integer tensor input.
            # Only the internal exact transfer accumulator changes representation.
            if environment and not isinstance(next(iter(environment.values())), gmpy2.mpz):
                environment = {k:gmpy2.mpz(v) for k,v in environment.items()}
            return original_step(self, environment, site, operator)
        modules[4].State.step = integer_step
    read = lambda p: json.loads(p.read_text())
    started = time.monotonic()
    expected = read(case/candidate/'exact/interval.json')
    lower = modules[2].check(read(case/'fixture.json'),
                           read(case/candidate/'exact/certificate.json'),
                           read(case/'nonsinglet.json'))
    lower_seconds = time.monotonic()-started
    upper = modules[4].check(read(case/'fixture.json'), read(case/'mps/state.json'))
    elapsed = time.monotonic()-started
    for got, field in [(lower['lower'], 'lower_Ha'), (upper['upper_Ha'], 'upper_Ha')]:
        if Fraction(got) != Fraction(expected[field]):
            raise AssertionError('Exact endpoint differs from stdlib checker')
    output = dict(gmpy2=gmpy2.version(), gmp_integer_transfers=integers, substituted_modules=substituted, lower_seconds=lower_seconds,
                  complete_replay_seconds=elapsed, exact_endpoints_match_stdlib=True,
                  authoritative_receipt='exact/interval.json', status='experimental_exact_arithmetic_backend')
    with (case/candidate/('gmp_integer_probe.json' if integers else 'gmp_probe.json')).open('x') as f:
        json.dump(output, f, indent=2)
    print(json.dumps(output), flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('case', type=Path)
    p.add_argument('candidate')
    p.add_argument('--integers', action='store_true')
    a = p.parse_args()
    run(a.case.resolve(), a.candidate, a.integers)
