"""Thin callers for existing state, nonsinglet, and exact accepting modules."""
import argparse
from fractions import Fraction as F
import json
from pathlib import Path
import time
from research.transfer_solver_20260915.budget import dump


def run(case, action):
    read = lambda path: json.loads(path.read_text())
    data = read(case/'fixture.json')
    if action == 'state':
        from research.correlated_pair_20260913.mps_spatial import run as discover
        out = case/'mps'
        out.mkdir(exist_ok=False)
        bond = min(144, 4**(data['modes']//4))
        discover(data, [data['particles']//2]*2, out, bond=bond, sweeps=10, initial=None)
    elif action == 'upper':
        from research.correlated_pair_20260913.mps_exact import check
        dump(case/'upper.json', check(data, read(case/'mps/state.json')))
    elif action == 'nonsinglet':
        from research.collective_completion_20260914.prepare import build
        from research.collective_completion_20260914 import search
        if min(data['particles'], data['modes']-data['particles']) < 2:
            raise ValueError('This campaign requires a nonempty nonsinglet screen')
        prepared = case/'nonsinglet_prepared'
        build(case.name, case/'fixture.json', case/'mps/state.json', prepared)
        search.OUT = case/'nonsinglet_search'
        search.solve(case.name, 4, seconds=90, spin=True, tag='_nonsinglet', prepared=prepared, mag=1)
        source = search.OUT/'candidates'/f'{case.name}_spin_r4_nonsinglet/round_0/certificate.json'
        dump(case/'nonsinglet.json', read(source))
    elif action == 'replay':
        import sys
        from research.collective_completion_20260914.spin_screen import check
        started = time.monotonic()
        cert = read(case/'stage2/export/certificate.json')
        removed, allowance = 0, F(0)
        blocks = []
        for block in cert['core']['blocks']:
            rows = []
            for row in block['factor']:
                total = sum(map(abs, row))
                if total < 1000:
                    removed += 1
                    allowance += F(total*total, cert['core']['denominator']**2)
                else:
                    rows.append(row)
            if rows:
                blocks.append({**block, 'factor': rows})
        cert['core']['blocks'] = blocks
        dump(case/'certificate.json', cert)
        receipt = check(data, cert, read(case/'nonsinglet.json'))
        dump(case/'lower.json', receipt)
        U, L = F(read(case/'upper.json')['upper_Ha']), F(receipt['lower'])
        if U < L:
            raise ValueError('Inconsistent exact endpoints')
        forbidden = [n for n in ('numpy', 'scipy', 'numba', 'quimb', 'pyscf', 'cvxpy') if n in sys.modules]
        if forbidden:
            raise AssertionError(('Numerical accepting dependency', forbidden))
        dump(case/'interval.json', {'lower_Ha': str(L), 'upper_Ha': str(U), 'width_Ha': str(U-L),
                                   'width_mHa': float(1000*(U-L)), 'target_met': U-L <= F(1, 625),
                                   'valid_on': receipt['valid_on'], 'removed_rows': removed,
                                   'removed_square_norm_bound_Ha': str(allowance),
                                   'lower_replay_seconds': time.monotonic()-started,
                                   'accepting_numerical_imports': forbidden})
        print(json.dumps({'width_mHa': float(1000*(U-L)), 'target_met': U-L <= F(1, 625)}), flush=True)
    else:
        raise ValueError('Unknown action')


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('case', type=Path)
    p.add_argument('action', choices=('state', 'upper', 'nonsinglet', 'replay'))
    a = p.parse_args()
    run(a.case.resolve(), a.action)
