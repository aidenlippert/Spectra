"""Accept an acceleration candidate using the unchanged, exact CPU checkers."""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import sys
import time


def replay(case, candidate):
    started = time.monotonic()
    out = case/candidate/'exact'
    out.mkdir(exist_ok=False)
    read = lambda p: json.loads(p.read_text())
    meta = read(case/'prepared/frame.json')
    for path, field in [(case/'fixture.json', 'fixture_sha256'),
                        (case/'mps/state.json', 'state_sha256')]:
        if hashlib.sha256(path.read_bytes()).hexdigest() != meta[field]:
            raise ValueError('Changed frozen input: '+path.name)
    from research.correlated_pair_20260913.mps_exact import check as upper_check
    from research.collective_completion_20260914.spin_screen import check as lower_check
    data = read(case/'fixture.json')
    upper = upper_check(data, read(case/'mps/state.json'))
    if F(upper['upper_Ha']) != F(read(case/'upper.json')['upper_Ha']):
        raise ValueError('Frozen upper did not replay')
    upper_seconds = time.monotonic()-started
    cert = read(case/candidate/'export/certificate.json')
    # Same row-pruning policy as the frozen transfer runner. The full residual
    # is reconstructed after pruning by the original accepting checker.
    removed, allowance, blocks = 0, F(0), []
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
    (out/'certificate.json').write_text(json.dumps(cert, separators=(',', ':'))+'\n')
    lower = lower_check(data, cert, read(case/'nonsinglet.json'))
    U, L = F(upper['upper_Ha']), F(lower['lower'])
    if U < L:
        raise ValueError('Inconsistent exact endpoints')
    forbidden = [n for n in ('numpy', 'scipy', 'cupy', 'torch', 'numba', 'quimb', 'pyscf', 'cvxpy') if n in sys.modules]
    if forbidden:
        raise AssertionError(('Numerical accepting dependency', forbidden))
    receipt = dict(lower_Ha=str(L), upper_Ha=str(U), width_Ha=str(U-L),
                   width_mHa=float(1000*(U-L)), target_met=U-L <= F(1,625),
                   removed_rows=removed, removed_square_norm_bound_Ha=str(allowance),
                   upper_replay_seconds=upper_seconds,
                   complete_replay_seconds=time.monotonic()-started,
                   accepting_numerical_imports=forbidden, valid_on=lower['valid_on'])
    for name, value in [('upper', upper), ('lower', lower), ('interval', receipt)]:
        (out/f'{name}.json').write_text(json.dumps(value, indent=2)+'\n')
    print(json.dumps({k:receipt[k] for k in ('width_mHa','target_met','complete_replay_seconds')}), flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('case', type=Path)
    p.add_argument('candidate')
    a = p.parse_args()
    replay(a.case.resolve(), a.candidate)
