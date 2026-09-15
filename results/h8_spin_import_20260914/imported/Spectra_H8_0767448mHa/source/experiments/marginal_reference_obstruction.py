"""Exact evidence that a selected determinant space cannot supply an accurate upper."""
from fractions import Fraction as F
import json
import math
from pathlib import Path
import time

from experiments.marginal_determinant_tree import DeterminantOracle
from experiments.marginal_schur_transfer import ldl_pivots
from experiments.marginal_implicit_certificate import rational_text

ROOT = Path(__file__).resolve().parents[1]


def replay(certificate):
    if certificate.get('kind') != 'retained_space_obstruction_v1':
        raise ValueError('Unsupported retained-space obstruction')
    oracle = DeterminantOracle(certificate)
    retained = certificate.get('retained_states')
    oracle.retained(retained)
    floor = F(certificate['retained_floor'])
    block, _, _ = oracle.retained_data(retained, F(0))
    shifted = [[x - floor * (i == j) for j, x in enumerate(row)] for i, row in enumerate(block)]
    pivots = ldl_pivots(shifted)
    if pivots is None:
        raise ValueError('Retained-space floor failed exact positivity')
    upper = oracle.upper(certificate.get('independent_upper'))
    if floor <= upper:
        raise ValueError('No positive retained-space obstruction established')
    return {'retained_floor': rational_text(floor), 'global_variational_upper': rational_text(upper),
            'retained_error_lower_bound': rational_text(floor - upper),
            'retained_error_lower_bound_float': float(floor - upper),
            'retained_dimension': len(retained), 'sector_dimension': math.comb(oracle.modes, oracle.particles),
            'upper_support_size': sum(bool(x) for x in certificate['independent_upper']['amplitudes']),
            'unique_action_states': len(oracle.cache), 'referenced_determinants': oracle.referenced_state_count(),
            'retained_pivots': [rational_text(x) for x in pivots],
            'scope': 'Every normalized state entirely in this retained space is at least the reported amount above the global ground energy of the rational Hamiltonian; no global lower endpoint is claimed'}


def run():
    root = ROOT / 'results/marginal_h6'
    out = root / 'retained_obstruction'
    if out.exists():
        raise ValueError('Preserve previous obstruction certificate')
    started = time.monotonic()
    fixture = json.loads((root / 'fixture.json').read_text())
    proposal = json.loads((root / 'sparse_p32/failure.json').read_text())
    reference = json.loads((root / 'reference_upper.json').read_text())
    oracle = DeterminantOracle(fixture)
    retained = proposal['retained_states']
    own_upper = oracle.upper(proposal['independent_upper'])
    block, _, _ = oracle.retained_data(retained, F(0))
    for margin in range(10):
        floor = F(math.floor((float(own_upper) - 10. ** (margin - 10)) * 10**12), 10**12)
        if ldl_pivots([[x - floor * (i == j) for j, x in enumerate(row)] for i, row in enumerate(block)]) is not None:
            break
    else:
        raise ValueError('No exact retained floor accepted')
    c = {key: fixture[key] for key in ('modes', 'particles', 'hamiltonian')}
    c.update(kind='retained_space_obstruction_v1', retained_states=retained,
             retained_floor=rational_text(floor), independent_upper=reference['independent_upper'])
    receipt = replay(c)
    receipt.update(elapsed_seconds=time.monotonic() - started,
                   selected_upper=rational_text(own_upper), source_proposal='results/marginal_h6/sparse_p32/failure.json')
    out.mkdir(parents=True)
    (out / 'certificate.json').write_text(json.dumps(c, indent=2) + '\n')
    (out / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({key: receipt[key] for key in ('retained_dimension', 'retained_error_lower_bound_float', 'elapsed_seconds')}), flush=True)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify')
    args = parser.parse_args()
    if args.verify:
        print(json.dumps(replay(json.loads(Path(args.verify).read_text())), indent=2))
    else:
        run()
