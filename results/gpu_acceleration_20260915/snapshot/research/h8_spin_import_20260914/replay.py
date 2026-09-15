"""Recompute the imported bound using only the original local accepting code."""
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import sys
import time

from research.correlated_pair_20260913.mps_exact import check as check_upper
from research.collective_completion_20260914.spin_screen import check as check_lower
from research.h8_spin_import_20260914.budget import ROOT, OUT


def dump(path, data):
    with path.open('x') as stream:
        json.dump(data, stream, indent=2)
        stream.write('\n')


def run():
    if (OUT/'manifest.json').exists():
        raise RuntimeError('Sealed validation pass')
    start = time.monotonic()
    destination = OUT/'local_exact_replay'
    destination.mkdir(exist_ok=False)
    package = OUT/'imported/Spectra_H8_0767448mHa'
    paths = {
        'fixture': ROOT/'results/certificate_scaling/active_space_ladder/h8/fixture.json',
        'state': ROOT/'results/correlated_pair_20260913/mps/h8_spatial_warm144/state.json',
        'singlet': package/'certificates/singlet.json',
        'nonsinglet': ROOT/'results/collective_completion_20260914/candidates/h8_spin_r1_nonsinglet/round_0/certificate.json',
    }
    read = lambda path: json.loads(path.read_text())
    data = {name: read(path) for name, path in paths.items()}
    frozen = read(ROOT/'results/sector_quotient_20260914/frozen_inputs.json')['h8']
    claim = read(package/'RESULT.json')
    print('Recomputing the original rational MPS expectation with local code.', flush=True)
    upper = check_upper(data['fixture'], data['state'])
    dump(destination/'upper.json', upper)
    print(json.dumps({'stage': 'upper_accepted', 'seconds': time.monotonic()-start,
                      'matches_frozen_upper': F(upper['upper_Ha']) == F(frozen['upper_Ha'])}), flush=True)
    lower = check_lower(data['fixture'], data['singlet'], data['nonsinglet'])
    dump(destination/'lower.json', lower)
    U, L = F(upper['upper_Ha']), F(lower['lower'])
    if U < L:
        raise ValueError('Inconsistent exact endpoints')
    if U != F(frozen['upper_Ha']):
        raise ValueError('Changed upper endpoint')
    if L != F(claim['lower_Ha']) or U != F(claim['upper_Ha']) or U-L != F(claim['width_Ha']):
        raise ValueError('Computed interval differs from the supplied claim')
    forbidden = [name for name in ('numpy', 'scipy', 'cvxpy', 'quimb', 'pyscf', 'numba', 'sympy') if name in sys.modules]
    if forbidden:
        raise AssertionError(('Numerical acceptance dependency', forbidden))
    modules = {}
    for name, module in tuple(sys.modules.items()):
        if not name.startswith(('research.', 'experiments.')):
            continue
        source = getattr(module, '__file__', None)
        if source is None:
            continue
        path = Path(source).resolve()
        if not path.is_relative_to(ROOT) or path.is_relative_to(OUT/'imported'):
            raise AssertionError(('Imported accepting code outside the original project', str(path)))
        modules[name] = {'path': str(path.relative_to(ROOT)), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}
    old = read(ROOT/'results/sector_quotient_20260914/candidates/linear_closure_repaired/interval.json')
    record = {'kind': 'independent_local_import_replay_v1', 'lower_Ha': str(L), 'upper_Ha': str(U),
              'width_Ha': str(U-L), 'width_mHa': float(1000*(U-L)), 'target_met': U-L <= F(1, 625),
              'matches_supplied_claim_exactly': True, 'matches_frozen_upper_exactly': True,
              'lower_improvement_vs_previous_compact_Ha': str(L-F(old['lower'])),
              'lower_improvement_vs_previous_strong_Ha': str(L-F(frozen['lower_Ha'])),
              'valid_on': lower['valid_on'], 'spin_defect_Ha': lower['singlet']['original_H_spin_defect_Ha'],
              'residual_wedge_witness_used': 'residual_wedge_witness' in data['singlet'],
              'upper_replay_seconds': upper['replay_seconds'], 'lower_replay_seconds': lower['replay_seconds'],
              'total_replay_seconds': time.monotonic()-start, 'numerical_libraries_imported': forbidden,
              'input_sha256': {name: {'path': str(path.relative_to(ROOT)), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}
                               for name, path in paths.items()},
              'local_accepting_modules': modules}
    dump(destination/'complete.json', record)
    print(json.dumps({k: v for k, v in record.items() if k not in ('input_sha256', 'local_accepting_modules', 'upper_Ha', 'width_Ha')}, indent=2))
    if not record['target_met']:
        raise SystemExit('The complete interval misses the target')


if __name__ == '__main__':
    run()
