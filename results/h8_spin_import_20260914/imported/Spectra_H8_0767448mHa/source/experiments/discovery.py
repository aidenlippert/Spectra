"""Finite Pauli support identification with a declared ideal derivative sensor.

The generic vocabulary is public. The observation oracle alone sees truth.
Physical derivative readout and its bounded error are assumptions, not apparatus.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction
from itertools import product
from time import perf_counter
import numpy as np
from experiments.pauli import commutator_i, multiply


@dataclass(frozen=True)
class Mechanism:
    support: tuple[str, ...]
    coefficients: tuple[Fraction, ...]

    def hamiltonian(self):
        return dict(zip(self.support, self.coefficients, strict=True))


def dictionary(n=3, max_weight=2):
    if not isinstance(n, int) or n < 1 or n > 6:
        raise ValueError('bounded dictionary requires 1 <= n <= 6')
    if not isinstance(max_weight, int) or not 1 <= max_weight <= n:
        raise ValueError('invalid maximum weight')
    return tuple(''.join(p) for p in product('IXYZ', repeat=n)
                 if 0 < sum(c != 'I' for c in p) <= max_weight)


def probe_bank(n=3, measurement_weight=1):
    """All Q of one declared weight, all nonidentity Pauli-biased states S."""
    observables = [p for p in dictionary(n, n)
                   if sum(c != 'I' for c in p) == measurement_weight]
    return tuple((q, s) for q in observables for s in dictionary(n, n))


def design(rows, candidates):
    """Tr((I+S)/2^n i[P,Q]); each row has at most one nonzero column."""
    a = np.zeros((len(rows), len(candidates)))
    index = {p: j for j, p in enumerate(candidates)}
    for i, (q, s) in enumerate(rows):
        phase, p = multiply(q, s)
        if phase in (1j, -1j) and p in index:
            a[i, index[p]] = int(commutator_i({p: 1}, q).get(s, 0))
    return a


def observe(truth, rows, rng, noise=Fraction(1, 500)):
    """Simulator boundary: exact rational signal plus bounded rational noise."""
    if noise < 0:
        raise ValueError('negative sensor noise')
    h = truth.hamiltonian()
    derivatives = {q: commutator_i(h, q) for q, _ in rows}
    values = []
    for q, s in rows:
        error = noise * Fraction(int(rng.integers(-1000, 1001)), 1000)
        values.append(derivatives[q].get(s, Fraction(0)) + error)
    return tuple(values)


def _fit(a, y, selected, ledger):
    ledger['least_squares_solves'] += 1
    ledger['least_squares_column_sum'] += len(selected)
    if not selected:
        return np.empty(0), y.copy()
    coef = np.linalg.lstsq(a[:, selected], y, rcond=None)[0]
    return coef, y - a[:, selected] @ coef


def fit_model(a, y, candidates, method, known=None, tolerance=0.006):
    """Only public design, observations and prior inferred coefficients enter."""
    known = known or {}
    ledger = {'correlation_dot_products': 0, 'correlation_multiply_adds': 0,
              'least_squares_solves': 0, 'least_squares_column_sum': 0}
    start = perf_counter()
    norms = np.sum(a * a, axis=0)
    if method == 'frozen':
        model = dict(known)
    elif method in ('retrieval', 'full_dictionary'):
        selected = (list(range(len(candidates))) if method == 'full_dictionary'
                    else [i for i, p in enumerate(candidates) if p in known])
        coef, _ = _fit(a, y, selected, ledger)
        model = {candidates[i]: float(c) for i, c in zip(selected, coef)}
    else:
        selected = ([i for i, p in enumerate(candidates) if p in known]
                    if method == 'sequential' else [])
        for _ in range(len(candidates) + 1):
            coef, residual = _fit(a, y, selected, ledger)
            if np.max(np.abs(residual), initial=0) <= tolerance:
                break
            scores = np.abs(a.T @ residual) / np.maximum(norms, 1)
            ledger['correlation_dot_products'] += len(candidates)
            ledger['correlation_multiply_adds'] += a.shape[0] * a.shape[1]
            scores[selected] = -1
            scores[norms == 0] = -1
            best = int(np.argmax(scores))
            if scores[best] <= tolerance / 2:
                break
            selected.append(best)
        coef, _ = _fit(a, y, selected, ledger)
        model = {candidates[i]: float(c) for i, c in zip(selected, coef)
                 if abs(c) > tolerance / 2}
    vector = np.array([model.get(p, 0.0) for p in candidates])
    residual = y - a @ vector
    identified = bool(np.all(norms > 0))  # independent one-column rows
    ledger['wall_seconds'] = perf_counter() - start
    return model, {'fit_rmse': float(np.sqrt(np.mean(residual ** 2))),
                   'fit_max_error': float(np.max(np.abs(residual), initial=0)),
                   'abstained': not identified or bool(np.max(np.abs(residual), initial=0) > tolerance),
                   'identifiable_dictionary': identified, 'cost': ledger}


def parameter_intervals(a, y_exact, noise):
    """Exact coverage conditional on the declared dictionary and sensor bound.

No confidence over unmodelled physics is asserted. Zero-feature rows can reject
some out-of-family processes; passing them cannot exclude all misspecification.
"""
    intervals = [None] * a.shape[1]
    for row, value in zip(a, y_exact, strict=True):
        nonzero = np.flatnonzero(row)
        if len(nonzero) > 1:
            raise ValueError('interval checker requires isolating Pauli probes')
        if not len(nonzero):
            if abs(value) > noise:
                return None, 'out_of_family_or_noise_bound_violated'
            continue
        j = int(nonzero[0]); scale = Fraction(int(row[j]))
        low, high = sorted(((value - noise) / scale, (value + noise) / scale))
        if intervals[j] is not None:
            low = max(low, intervals[j][0]); high = min(high, intervals[j][1])
        if low > high:
            return None, 'inconsistent_sensor_intervals'
        intervals[j] = (low, high)
    if any(v is None for v in intervals):
        return None, 'unidentified_coefficients'
    return intervals, 'bounded_within_declared_family'


def _truths(rng, candidates):
    # This is the simulator's generator; fit_model never receives its labels.
    base = ('XII', 'ZII')
    first_pool = [p for p in candidates if p[0] != 'I' and p[1] != 'I' and p[2] == 'I'
                  and any(commutator_i({p: 1}, q) for q in base)]
    one = first_pool[int(rng.integers(len(first_pool)))]
    second_pool = [p for p in candidates if p[2] != 'I' and sum(c != 'I' for c in p) == 2
                   and commutator_i({p: 1}, one)]
    two = second_pool[int(rng.integers(len(second_pool)))]
    return [Mechanism(base, (Fraction(7, 10), Fraction(-9, 20))),
            Mechanism(base + (one,), (Fraction(7, 10), Fraction(-9, 20), Fraction(11, 20))),
            Mechanism(base + (one, two), (Fraction(7, 10), Fraction(-9, 20), Fraction(11, 20), Fraction(-2, 5)))]


def benchmark(seed=7, noise=Fraction(1, 500)):
    rng = np.random.default_rng(seed)
    candidates = list(dictionary())
    rng.shuffle(candidates)
    candidates = tuple(candidates)
    rows = list(probe_bank(measurement_weight=1)); rng.shuffle(rows); rows = tuple(rows)
    test_rows = probe_bank(measurement_weight=2)  # disjoint Q, actual new probes
    a, at = design(rows, candidates), design(test_rows, candidates)
    truths = _truths(rng, candidates)
    prev = {}; frozen = None; stages = []; cumulative = {}
    methods = ('frozen', 'retrieval', 'sequential', 'from_scratch', 'full_dictionary')
    for k, truth in enumerate(truths):
        # Each stage has a genuinely new Hamiltonian context for reuse testing.
        contexts = [truth, Mechanism(truth.support, tuple(c * Fraction(13 + i + k, 10)
                                                         for i, c in enumerate(truth.coefficients)))]
        context_results = []
        for context_index, current in enumerate(contexts):
            exact = observe(current, rows, rng, noise)
            test = observe(current, test_rows, rng, noise)
            y, yt = np.array([float(v) for v in exact]), np.array([float(v) for v in test])
            intervals, status = parameter_intervals(a, exact, noise)
            measured = {}
            # Before the first discovery all adaptive methods start empty.
            for method in methods:
                effective = ('from_scratch' if k == 0 and context_index == 0 and method in ('frozen', 'retrieval') else method)
                old = frozen if method == 'frozen' and frozen is not None else prev
                model, metrics = fit_model(a, y, candidates, effective, old)
                vector = np.array([model.get(p, 0) for p in candidates])
                err = yt - at @ vector
                metrics.update({'holdout_rmse': float(np.sqrt(np.mean(err ** 2))),
                                'holdout_max_error': float(np.max(np.abs(err))),
                                'support': sorted(p for p, c in model.items() if abs(c) > 0.003),
                                'coefficients': model,
                                'oracle_support_correct': set(p for p, c in model.items() if abs(c) > 0.003) == set(current.support),
                                'sensor_records_used': len(rows),
                                'validation_records': len(test_rows),
                                'public_design_entries': a.size,
                                'public_feature_builds': len(rows)})
                if metrics['holdout_max_error'] > 0.006:
                    metrics['abstained'] = True
                # Costs shared by all methods are shown per method, not free.
                subtotal = metrics['cost']['correlation_multiply_adds']
                cumulative[method] = cumulative.get(method, 0) + subtotal
                metrics['cost']['cumulative_correlation_multiply_adds'] = cumulative[method]
                measured[method] = metrics
            seq = measured['sequential']['coefficients']
            rational = {p: Fraction(v).limit_denominator(10000) for p, v in seq.items()}
            if intervals is None:
                eta = None
            else:
                eta = sum(max(abs(rational.get(p, 0) - lo), abs(rational.get(p, 0) - hi))
                          for p, (lo, hi) in zip(candidates, intervals))
            context_results.append({'kind': 'discovery' if context_index == 0 else 'new_coefficients',
                'methods': measured, 'interval_status': status,
                'conditional_model_norm_error': None if eta is None else str(eta),
                'learned_hamiltonian': {p: str(v) for p, v in rational.items()},
                'oracle_hamiltonian': {p: str(c) for p, c in current.hamiltonian().items()},
                'parameter_intervals': None if intervals is None else {p: [str(lo), str(hi)] for p, (lo, hi) in zip(candidates, intervals)}})
            if context_index == 0:
                prev = dict(seq)
                if frozen is None:
                    frozen = dict(seq)
        stages.append({'stage': k, 'contexts': context_results})
    # Out-of-class third-order interaction: diagnostic rows must reject it.
    outside = Mechanism(('XXX',), (Fraction(1, 2),))
    outside_y = observe(outside, rows, rng, noise)
    _, outside_status = parameter_intervals(a, outside_y, noise)
    zero_a = np.zeros((4, len(candidates)))
    _, unident = fit_model(zero_a, np.zeros(4), candidates, 'from_scratch')
    return {'seed': seed, 'n_qubits': 3, 'dictionary': list(candidates),
            'noise_absolute_bound': str(noise), 'sensor': 'ideal_initial_time_pauli_derivative',
            'sensor_definition': 'rho=(I+S)/2^n; y=Tr(rho i[H,Q]); bounded rational additive error',
            'identity_energy_offset': 'fixed zero; not observable from commutators',
            'stages': stages, 'misspecified_status': outside_status,
            'unidentifiable_abstained': unident['abstained'],
            'claim': 'closed-vocabulary support recovery; no mechanism vocabulary or experimental strategy learning',
            'resource_scope': 'counts derivative-oracle records, algebra and fits; physical shots/derivative-estimation time are not modelled',
            'compounding_demonstrated': False}
