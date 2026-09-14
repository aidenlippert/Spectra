"""Evaluate the proved full-enrichment error bound with rational inputs."""
from fractions import Fraction as F
from pathlib import Path
import json
import math

from experiments.marginal_schur_transfer import complement_lower
from experiments.marginal_general_schur import fixture, norm_bound
from experiments.marginal_collective import hopping_polynomial
from experiments.marginal_symbolic import add, scale

root = Path(__file__).resolve().parent
out = root / 'enrichment_order_bound.json'
if out.exists():
    raise ValueError('Preserve previous bound calculation')
rows = []
for name, mixed, folder in [('cycle', False, 'cycle_1_100_targeted_upper2_bounded'),
                            ('mixed', True, 'mixed_1_100_targeted_grouped_upper1')]:
    receipt = json.loads((root / folder / 'independent_replay.json').read_text())
    exact_u, exact_b = F(receipt['upper']), F(receipt['lower'])
    c0 = complement_lower()
    denominator = 10**9
    c = F(math.floor(c0 * denominator), denominator)
    u = F(math.ceil(exact_u * denominator), denominator)
    b = F(math.floor(exact_b * denominator), denominator)
    assert c <= c0 and b <= exact_b and u >= exact_u
    delta = add(fixture(F(1, 100), mixed), scale(hopping_polynomial(10, F(1, 5)), -1))
    eta = norm_bound(delta, 'hermitian_pairs')
    spectral_upper = 11 + eta
    rho, width = eta / (c - u), spectral_upper - b
    assert 0 <= rho < 1
    target = F(1, 10**7)
    for order in range(100):
        error = width * rho ** (2 * order + 2)
        if error <= target:
            break
    else:
        raise ValueError('Order search budget exhausted')
    rows.append({'fixture': name, 'source_replay': str((root / folder / 'independent_replay.json').relative_to(root.parent.parent)),
                 'complement_lower': str(c), 'perturbation_norm_bound': str(eta),
                 'energy_lower': str(b), 'energy_upper': str(u), 'spectral_upper': str(spectral_upper),
                 'rho': str(rho), 'rho_float': float(rho), 'target_error': str(target),
                 'sufficient_full_enrichment_order': order, 'guaranteed_ritz_error': str(error),
                 'guaranteed_ritz_error_float': float(error),
                 'scope': 'Full H0-closed enrichment chain only; not a guarantee for selected targeting directions'})
out.write_text(json.dumps(rows, indent=2) + '\n')
for row in rows:
    print(row['fixture'], row['rho_float'], row['sufficient_full_enrichment_order'], row['guaranteed_ritz_error_float'])
