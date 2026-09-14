"""Finite representation revision with sampled-only access and mismatch checks.

This is a separate experiment from the complementary-risk theorem. It does not
claim a warm-start sample advantage: changing systems are fully remeasured.
"""
from dataclasses import asdict, dataclass
from fractions import Fraction as F
from itertools import combinations, product
from pathlib import Path
import json
import numpy as np
from experiments.v2_laws import CalibrationRecord, learn, required_shots, verify_certificate, _exact_union_bound


def _bit_tuple(x, d):
    if type(x) is not tuple or len(x) != d or any(type(v) is not int or v not in (0, 1) for v in x):
        raise ValueError('strict binary context required')


@dataclass(frozen=True)
class Law:
    linear: tuple
    quadratic: tuple = ()
    cubic: tuple = ()
    constant: int = 0

    def value(self, x):
        d = len(self.linear)
        _bit_tuple(self.linear, d)
        _bit_tuple(x, d)
        if type(self.constant) is not int or self.constant not in (0, 1):
            raise ValueError('binary constant required')
        value = (self.constant + sum(a * b for a, b in zip(x, self.linear))) % 2
        for term in self.quadratic + self.cubic:
            if len(term) not in (2, 3) or tuple(sorted(set(term))) != term or any(type(i) is not int or not 0 <= i < d for i in term):
                raise ValueError('invalid monomial')
            value ^= int(all(x[i] for i in term))
        return value


def features(x):
    return (1,) + x + tuple(x[i] & x[j] for i, j in combinations(range(len(x)), 2))


def sample_callback(law, eta, rng):
    if not isinstance(eta, F) or not 0 <= eta < F(1, 2):
        raise ValueError('exact crossover below 1/2 required')
    def callback(x, n):
        truth = law.value(x)
        return tuple(truth ^ int(t < eta.numerator) for t in rng.integers(0, eta.denominator, size=n))
    return callback


def acquire(callback, d, eta=F(1, 10), delta=F(1, 1000), shots=None):
    if type(d) is not int or not 3 <= d <= 8:
        raise ValueError('d must be in 3..8')
    training = tuple(x for x in product((0, 1), repeat=d) if sum(x) <= 2)
    n = required_shots(len(training), eta, delta) if shots is None else shots
    if type(n) is not int or n < 1 or n % 2 == 0:
        raise ValueError('positive odd shot count required')
    if _exact_union_bound(len(training), n, eta) > delta:
        return {'status': 'insufficient_shots', 'source_shots': 0, 'model': None}
    # Query zero and units first, then pair interventions. Callback returns only
    # noisy bits; coefficients and source implementation are inaccessible here.
    ordered = tuple(sorted(training, key=lambda x: (sum(x), x)))
    records = tuple(CalibrationRecord(x, tuple(callback(x, n))) for x in ordered)
    labels = {rec.context: int(sum(rec.outcomes) > n // 2) for rec in records}
    zero = (0,) * d
    constant = labels[zero]
    units = tuple(tuple(int(i == j) for i in range(d)) for j in range(d))
    linear = tuple(labels[x] ^ constant for x in units)
    residuals = []
    support = []
    for i, j in combinations(range(d), 2):
        x = tuple(int(k in (i, j)) for k in range(d))
        bit = labels[x] ^ labels[units[i]] ^ labels[units[j]] ^ constant
        residuals.append({'pair': (i, j), 'context': x, 'residual': bit})
        if bit:
            support.append((i, j))
    # Independent full-feature conventional linear solve over the identical data.
    expanded = tuple(CalibrationRecord(features(r.context), r.outcomes) for r in records)
    conventional = learn(expanded, eta, delta)
    if not verify_certificate(expanded, conventional, eta, delta):
        raise AssertionError('full-feature certificate failed')
    fitted = Law(linear, tuple(support), constant=constant)
    expected_coefficients = (constant,) + linear + tuple(int(pair in support) for pair in combinations(range(d), 2))
    if conventional.mask != expected_coefficients:
        raise AssertionError('independent coefficient extraction disagrees')
    return {'status': 'model_conditional_certificate', 'model': fitted,
            'initial_representation': 'affine linear',
            'representation_after_pair_probes': 'quadratic' if support else 'affine linear',
            'residual_witnesses': [r for r in residuals if r['residual']],
            'records': records, 'conventional_certificate': conventional,
            'conventional_coefficients_tie': True,
            'exact_error_bound': str(_exact_union_bound(len(training), n, eta)),
            'contexts_queried': len(training), 'shots_per_context': n,
            'source_shots': len(training) * n,
            'degree_membership_established_by_calibration': False}


def diagnose(callback, acquired, d, eta=F(1, 10), delta=F(1, 1000)):
    if acquired['model'] is None:
        return {'status': 'abstain', 'source_shots': 0}
    heldout = tuple(x for x in product((0, 1), repeat=d) if sum(x) >= 3)
    n = required_shots(len(heldout), eta, delta)
    records = tuple(CalibrationRecord(x, tuple(callback(x, n))) for x in heldout)
    disagreements = [r.context for r in records
                     if int(sum(r.outcomes) > n // 2) != acquired['model'].value(r.context)]
    return {'status': 'model_mismatch_abstain' if disagreements else 'no_heldout_disagreement',
            'disagreements': disagreements, 'records': records,
            'contexts_checked': len(heldout), 'source_shots': n * len(heldout),
            'exact_label_error_bound': str(_exact_union_bound(len(heldout), n, eta)),
            'labels_used_for_fitting': False,
            'scope': 'all previously unqueried points of this finite binary cube; no guarantee outside that domain'}


def _serialize(acq):
    return {k: (asdict(v) if k in ('model', 'conventional_certificate') and v is not None
                else [asdict(r) for r in v] if k == 'records' else v) for k, v in acq.items()}


def run(seed=61, d=4):
    if type(d) is not int or not 3 <= d <= 8:
        raise ValueError('d must be in 3..8')
    rng = np.random.default_rng(seed)
    linear = tuple(int(b) for b in rng.integers(0, 2, size=d))
    laws = (Law(linear, ((0, 1),)), Law(linear, ((0, 1), (1, 2))))
    stages = []
    prior_support = set()
    for law in laws:
        callback = sample_callback(law, F(1, 10), rng)
        acq = acquire(callback, d)
        diagnostic = diagnose(callback, acq, d)
        recovered = acq['model'] == law
        new_terms = sorted(set(acq['model'].quadratic) - prior_support)
        prior_support = set(acq['model'].quadratic)
        stages.append({'acquisition': _serialize(acq), 'diagnostic': _serialize(diagnostic),
                       'evaluator_exact_law_recovery': recovered, 'new_interactions': new_terms})
    cubic = Law(linear, ((0, 1),), ((0, 1, 2),))
    cb = sample_callback(cubic, F(1, 10), rng)
    ca = acquire(cb, d)
    cd = diagnose(cb, ca, d)
    training = tuple(r.context for r in ca['records'])
    # Algebraic indistinguishability on training is checked independently of noise.
    indistinguishable = all(cubic.value(x) == laws[0].value(x) for x in training)
    result = {'protocol': {'seed': seed, 'dimension': d, 'noise': '1/10',
                            'delta_per_fit_or_diagnostic': '1/1000',
                            'training_weight': '<=2', 'diagnostic_weight': '>=3'},
              'stages': stages,
              'cubic_negative_control': {'training_indistinguishable_from_quadratic': indistinguishable,
                                         'acquisition': _serialize(ca), 'diagnostic': _serialize(cd)},
              'costs': {'warm_source_shots': sum(s['acquisition']['source_shots'] for s in stages),
                        'reset_source_shots': sum(s['acquisition']['source_shots'] for s in stages),
                        'warm_shot_reduction': 0,
                        'heldout_diagnostic_shots': sum(s['diagnostic']['source_shots'] for s in stages),
                        'negative_control_shots': ca['source_shots'] + cd['source_shots']},
              'claim': 'sampled interaction discovery, heldout transfer, mismatch detection; no measured sample compounding',
              'physical_experiments_performed': False}
    return result


if __name__ == '__main__':
    result = run()
    path = Path(__file__).resolve().parents[1] / 'results/v2/revision_results.json'
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'claim': result['claim'], 'costs': result['costs'],
                      'recovered': [s['evaluator_exact_law_recovery'] for s in result['stages']],
                      'cubic_status': result['cubic_negative_control']['diagnostic']['status']}, indent=2))
