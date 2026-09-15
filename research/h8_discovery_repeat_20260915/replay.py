"""Accept a newly discovered certificate with the original exact local checker."""
from fractions import Fraction as F
import json
import sys
import time

from research.correlated_pair_20260913.mps_exact import check as check_upper
from research.collective_completion_20260914.spin_screen import check as check_lower
from research.h8_discovery_repeat_20260915.run import ROOT, OUT, INPUTS, dump, sha


def run():
    started = time.monotonic()
    folder = OUT/'exact_replay'
    folder.mkdir(exist_ok=False)
    read = lambda p: json.loads(p.read_text())
    fixture = read(INPUTS['fixture'])
    print('Recomputing the exact rational MPS upper.', flush=True)
    upper = check_upper(fixture, read(INPUTS['state']))
    dump(folder/'upper.json', upper)
    print('Expanding the new singlet factors and checking both spin sectors.', flush=True)
    lower = check_lower(fixture, read(OUT/'certificate.json'), read(INPUTS['nonsinglet']))
    dump(folder/'lower.json', lower)
    U, L = F(upper['upper_Ha']), F(lower['lower'])
    frozen = read(INPUTS['frozen_endpoints'])['h8']
    if U < L or U != F(frozen['upper_Ha']):
        raise ValueError('Inconsistent endpoints or changed upper')
    forbidden = [n for n in ('numpy', 'scipy', 'cvxpy', 'quimb', 'pyscf', 'numba', 'sympy') if n in sys.modules]
    if forbidden:
        raise ValueError(('Numerical imports on the accepting path', forbidden))
    modules = {}
    for name, mod in tuple(sys.modules.items()):
        if name.startswith(('research.', 'experiments.')) and getattr(mod, '__file__', None):
            modules[name] = {'path': mod.__file__, 'sha256': sha(mod.__file__)}
    record = {'kind': 'fresh_H8_optimization_exact_replay_v1', 'lower_Ha': str(L), 'upper_Ha': str(U),
              'width_Ha': str(U-L), 'width_mHa': float(1000*(U-L)), 'target_met': U-L <= F(1, 625),
              'valid_on': lower['valid_on'], 'matches_frozen_upper': True,
              'original_H_spin_defect_Ha': lower['singlet']['original_H_spin_defect_Ha'],
              'certificate_sha256': sha(OUT/'certificate.json'),
              'factor_rows': lower['singlet']['factor_rows'], 'factor_nonzeros': lower['singlet']['factor_nonzeros'],
              'residual_l1_Ha': lower['singlet']['residual_l1'],
              'numerical_libraries_imported': forbidden, 'accepting_modules': modules,
              'upper_replay_seconds': upper['replay_seconds'], 'lower_replay_seconds': lower['replay_seconds'],
              'total_replay_seconds': time.monotonic()-started,
              'cold_discovery_claimed': False}
    dump(folder/'complete.json', record)
    print(json.dumps({k: v for k, v in record.items() if k not in ('upper_Ha', 'width_Ha', 'accepting_modules')}, indent=2))


if __name__ == '__main__':
    run()
