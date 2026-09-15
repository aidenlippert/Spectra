"""Common-target complementarity certificate and controlled retention study.

The mathematical certificate bounds optimal one-scalar Bayes risks; the
simulation evaluates explicit nulling policies and is only a cross-check.
"""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import json
import numpy as np
from experiments.v5_frame import C, S, Q, BASE_LAMBDAS, frame, null_probe, dot
from experiments.v5_core import (Contract, LearnedPrefix, acquire_next, make_source,
    ACTUATOR_BOUND, PREPARATION_TV, PREPARATION_TIME)

ROOT = Path(__file__).resolve().parents[1]
MASKS = ('neither', 'A', 'B', 'both')


def certificate():
    eta = ACTUATOR_BOUND
    probe = null_probe((1, 1))
    cases = []
    for prefix in product((-1, 1), repeat=2):
        modes, residual = frame(prefix)
        nuisance_lower = sum(F(lam) * max(abs(dot(probe, mode)) - eta, F(0))**2
                             for lam, mode in zip(BASE_LAMBDAS, modes))
        half_difference_upper = 2 * BASE_LAMBDAS[2] * C * S * (
            abs(dot(probe, residual)) + eta) * (abs(dot(probe, Q[3])) + eta)
        if prefix != (1, 1):
            assert nuisance_lower >= 125 * half_difference_upper
        cases.append({'prefix': prefix, 'nuisance_lower_per_alpha': str(nuisance_lower),
                      'half_difference_upper_per_alpha': str(half_difference_upper),
                      'mismatch_ratio_at_most_1_over_125': prefix != (1, 1)})
    # Spectral certificate for B-only side information: summing over a, the
    # residual Gram has eigenvalues 2*s^4, 2*(1-s^4), 0 in the nuisance space.
    residual_eigenvalues = (2*S**4, 2*(1-S**4), F(0))
    assert max(residual_eigenvalues) <= F(32, 25)
    assert all(2 - (1 + F(9, 16)) * ev >= 0 for ev in residual_eigenvalues)
    p = PREPARATION_TV
    single_lower = F(123, 496) - p
    both_upper = F(11, 500) + p
    neither_upper = F(1, 4) * F(11, 500) + F(3, 4) * (F(1, 2) + F(1, 248)) + p
    ideal_margin = 2 * single_lower - both_upper - neither_upper
    # A record can be wrong with probability e; B and (A,B) with 2e because
    # acquisition of B uses the possibly mistaken A observation.
    acquisition_correction = 5 * (F(1, 399) + p)
    acquired_margin = ideal_margin - acquisition_correction
    assert acquired_margin == F(2243, 24800) - F(5, 399) - 9*p
    assert acquired_margin > F(77, 1000)
    return {'side_information': 'only retained sign bits; no discarded probes/mode vectors',
            'target': 'one fresh independent third-mode sign, common across all four risks',
            'optimal_risk_A_lower': str(single_lower), 'optimal_risk_B_lower': str(single_lower),
            'optimal_risk_both_upper': str(both_upper), 'optimal_risk_neither_upper': str(neither_upper),
            'true_sign_complementarity_lower': str(ideal_margin),
            'acquisition_coupling_correction': str(acquisition_correction),
            'acquired_complementarity_lower': str(acquired_margin),
            'acquired_complementarity_lower_decimal': float(acquired_margin),
            'B_only_residual_gram_eigenvalues': [str(x) for x in residual_eigenvalues],
            'mismatch_certificates': cases}


def retained_policy(read, contract, mask, bits):
    """The only calibration input is the explicitly retained bit tuple."""
    if mask not in MASKS:
        raise ValueError('unknown retention mask')
    expected = {'neither': 0, 'A': 1, 'B': 1, 'both': 2}[mask]
    if type(bits) is not tuple or len(bits) != expected or any(type(x) is not int or x not in (-1, 1) for x in bits):
        raise ValueError('only permitted sign bits may be retained')
    contract.check()
    if len(contract.lineage) != 3:
        raise ValueError('common final-stage target required')
    guess = {'neither': (1, 1), 'A': (bits[0], 1) if bits else (1, 1),
             'B': (1, bits[0]) if bits else (1, 1),
             'both': bits if len(bits) == 2 else (1, 1)}[mask]
    return 1 - 2 * read(3, null_probe(guess))


def run(worlds=128, targets=16):
    if not 1 <= worlds <= 128 or not 1 <= targets <= 16:
        raise ValueError('bounded study exceeded')
    rows = []; losses = {mask: 0 for mask in MASKS}
    calibration_errors = 0
    for index in range(worlds):
        seed = 91001 + index
        rng = np.random.default_rng(seed)
        a, b = (int(x) for x in rng.choice((-1, 1), size=2))
        lineage = (f'cal-{index}-a', f'cal-{index}-b')
        # Calibration is completed before any target sign or target seed is drawn.
        contract = Contract(lineage, F(1))
        read, ledger, records = make_source((a, b), contract, seed * 101)
        first = acquire_next(read, Contract(lineage[:1], F(1)))
        second = acquire_next(read, contract, first)
        ah, bh = second.signs
        calibration_errors += int(second.signs != (a, b))
        row = {'seed': seed, 'evaluator_prefix': (a, b), 'acquired_bits': (ah, bh),
               'calibration_ledger': ledger, 'calibration_records': records, 'targets': []}
        alphas = rng.choice(np.arange(1001, 2000), size=targets, replace=False)
        for j in range(targets):
            c = int(rng.choice((-1, 1)))
            target_seed = int(rng.integers(1, 10**8))
            target_contract = Contract(lineage + (f'target-{index}-{j}',), F(int(alphas[j]), 1000))
            item = {'alpha': str(target_contract.alpha), 'evaluator_target': c,
                    'estimates': {}, 'records': {}, 'ledgers': {}}
            side_info = {'neither': (), 'A': (ah,), 'B': (bh,), 'both': (ah, bh)}
            for mask in MASKS:
                reader, target_ledger, target_records = make_source((a,b,c), target_contract, target_seed)
                estimate = retained_policy(reader, target_contract, mask, side_info[mask])
                item['estimates'][mask] = estimate
                item['records'][mask] = target_records
                item['ledgers'][mask] = target_ledger
                losses[mask] += int(estimate != c)
            row['targets'].append(item)
        rows.append(row)
    n = worlds * targets
    contrast = F(losses['A'] + losses['B'] - losses['both'] - losses['neither'], n)
    return {'certificate': certificate(), 'physical_experiment': False,
            'worlds': worlds, 'targets_per_world': targets, 'calibration_prefix_errors': calibration_errors,
            'acquisition_reads': 2*worlds, 'counterfactual_target_reads': 4*n,
            'total_relaxation_time': (2*worlds+4*n)*PREPARATION_TIME,
            'loss_counts': losses, 'empirical_policy_risks': {m: losses[m]/n for m in MASKS},
            'empirical_policy_complementarity': str(contrast),
            'interpretation': 'descriptive paired policy contrast, not an optimal-risk estimate or iid confidence interval',
            'rows': rows}


if __name__ == '__main__':
    output = run()
    path = ROOT / 'results/v5/complementarity_results.json'
    path.write_text(json.dumps(output, indent=2) + '\n')
    print(json.dumps({k:v for k,v in output.items() if k != 'rows'}, indent=2))
