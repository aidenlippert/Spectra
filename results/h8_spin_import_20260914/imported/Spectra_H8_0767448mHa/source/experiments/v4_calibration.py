"""Learn signed experiment moments from charged labelled calibration episodes.

The learner receives a sample callback, action count and tolerances. It never
receives a likelihood table, selected hidden model, or teacher action value.
The source must keep the hidden model fixed within each episode, including
repeated actions on conditionally independent copies.
"""
from fractions import Fraction as F
from itertools import product
from math import ceil
from pathlib import Path
import json
import numpy as np
from experiments.v4_moments import plan_from_moments, signed_gram_plan

ROOT = Path(__file__).resolve().parents[1]


def sample_budget(action_count, epsilon, delta):
    if type(action_count) is not int or not 1 <= action_count <= 12:
        raise ValueError('bounded action count required')
    if type(epsilon) is not F or not 0 < epsilon < 1 or type(delta) is not F or not 0 < delta < 1:
        raise ValueError('exact tolerances in (0,1) required')
    moments = 1 + action_count + action_count * (action_count + 1) // 2
    power = 0
    while 2 * moments * F(1, 2 ** power) > delta:
        power += 1
    # ln(2) < 1 gives exp(-power) <= 2^-power. This conservative
    # integer bound avoids making floating logarithm rounding a proof step.
    shots = ceil(F(2 * power) / epsilon ** 2)
    if moments * shots > 5_000_000:
        raise ValueError('calibration episode budget exceeds cap')
    return moments, shots, power


def learn_from_calibration(draw, action_count, epsilon=F(1, 40), delta=F(1, 100)):
    count, per_moment, power = sample_budget(action_count, epsilon, delta)
    work = {'labelled_episodes': 0, 'target_label_readouts': 0,
            'visible_readouts': 0, 'conditioned_visible_copy_preparations': 0,
            'signed_statistic_products': 0, 'moment_sums': 0}
    def estimate(actions):
        labels, outcomes = draw(tuple(actions), per_moment)
        labels, outcomes = np.asarray(labels), np.asarray(outcomes)
        if labels.shape != (per_moment,) or outcomes.shape != (per_moment, len(actions)):
            raise ValueError('malformed calibration record')
        if not np.isin(labels, (-1, 1)).all() or not np.isin(outcomes, (0, 1)).all():
            raise ValueError('invalid calibration observations')
        statistic = labels.astype(np.int64)
        for column in range(len(actions)):
            statistic = statistic * outcomes[:, column]
        work['labelled_episodes'] += per_moment
        work['target_label_readouts'] += per_moment
        work['visible_readouts'] += per_moment * len(actions)
        work['conditioned_visible_copy_preparations'] += per_moment * len(actions)
        work['signed_statistic_products'] += per_moment * len(actions)
        work['moment_sums'] += 1
        return F(int(statistic.sum()), per_moment)
    bias = estimate(())
    first = tuple(estimate((a,)) for a in range(action_count))
    gram = [[F(0) for _ in range(action_count)] for _ in range(action_count)]
    for a in range(action_count):
        for b in range(a, action_count):
            gram[a][b] = gram[b][a] = estimate((a, b))
    policy = policy_from_estimated_moments(bias, first, tuple(tuple(row) for row in gram), epsilon)
    return dict(policy, epsilon=epsilon, delta=delta,
                moment_count=count, episodes_per_moment=per_moment,
                union_failure_upper_bound=2 * count * F(1, 2 ** power), work=work,
                scope='conditional on labelled source, independent episodes and same hidden model per episode')


def policy_from_estimated_moments(bias, first, gram, epsilon):
    if type(epsilon) is not F or not 0 < epsilon < 1:
        raise ValueError('exact moment error radius required')
    plan = plan_from_moments(bias, first, gram)
    a = plan['first_action']; children = plan['child_actions'][a]
    decisions = []
    for y, b in enumerate(children):
        masses = ((bias - first[a] - first[b] + gram[a][b], first[b] - gram[a][b])
                  if y == 0 else (first[a] - gram[a][b], gram[a][b]))
        decisions.append(tuple(1 if value >= 0 else -1 for value in masses))
    kappa = F(9, 2) * epsilon
    surrogate_risk = plan['values'][a]
    lower, upper = max(F(0), surrogate_risk - kappa), min(F(1), surrogate_risk + kappa)
    return {'first_action': a, 'child_actions': children, 'leaf_decisions': tuple(decisions),
            'moments': plan['moments'], 'surrogate_risk': surrogate_risk,
            'accepted': lower <= upper,
            'conditional_risk_interval': (lower, upper) if lower <= upper else None,
            'conditional_regret_bound': 2 * kappa,
            'refusal_reason': None if lower <= upper else 'moment intervals incompatible with any risk in [0,1]'}


def labelled_source(table, signs, seed, *, redraw_each_action=False):
    """Evaluator-owned finite uniform-mixture source. Exposes only a callback."""
    m, actions = len(signs), len(table[0])
    denominator = 1
    from math import lcm
    for row in table:
        for p in row:
            denominator = lcm(denominator, p.denominator)
    numerators = np.asarray([[int(p * denominator) for p in row] for row in table], dtype=np.int64)
    labels = np.asarray(signs, dtype=np.int8)
    rng = np.random.default_rng(seed)
    def draw(selected, n):
        if not 0 <= len(selected) <= 2 or any(type(a) is not int or not 0 <= a < actions for a in selected):
            raise ValueError('invalid calibration action')
        hidden = rng.integers(0, m, size=n)
        outputs = np.empty((n, len(selected)), dtype=np.int8)
        for column, a in enumerate(selected):
            used = rng.integers(0, m, size=n) if redraw_each_action else hidden
            outputs[:, column] = (rng.integers(0, denominator, size=n) < numerators[used, a])
        return labels[hidden], outputs
    return draw


def exact_policy_risk(policy, prior, table, signs):
    """Evaluator-only direct joint-outcome integration; not used in learning."""
    risk = F(0); a = policy['first_action']
    for y, b in enumerate(policy['child_actions']):
        for z, decision in enumerate(policy['leaf_decisions'][y]):
            for i, p in enumerate(prior):
                likelihood = (table[i][a] if y else 1 - table[i][a]) * (table[i][b] if z else 1 - table[i][b])
                if decision != signs[i]: risk += p * likelihood
    return risk


def test_world(right=False, visibility=F(4, 5)):
    bits = tuple(product((0, 1), repeat=3))
    signs = tuple((-1) ** (u + v) for u, v, w in bits)
    axes = (2, 2, 0, 1) if right else (0, 1, 2, 2)
    table = tuple(tuple((1 + visibility * (2 * row[a] - 1)) / 2 for a in axes) for row in bits)
    return (F(1, 8),) * 8, table, signs


def serialize(value):
    if isinstance(value, F): return str(value)
    if isinstance(value, dict): return {k: serialize(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)): return [serialize(v) for v in value]
    return value


def run():
    rows = []
    for index in range(4):
        prior, table, signs = test_world(bool(index % 2))
        learned = learn_from_calibration(labelled_source(table, signs, 14001 + index), 4)
        risk = exact_policy_risk(learned, prior, table, signs)
        reference = signed_gram_plan(prior, table, signs)
        lower, upper = learned['conditional_risk_interval']
        actual_error = max(abs(learned['moments']['bias'] - reference['moments']['bias']),
                           max(abs(a-b) for a,b in zip(learned['moments']['first'], reference['moments']['first'])),
                           max(abs(learned['moments']['gram'][a][b] - reference['moments']['gram'][a][b]) for a in range(4) for b in range(4)))
        rows.append({'seed': 14001 + index, 'right_world': bool(index % 2),
                     'learned': learned, 'actual_risk': risk, 'optimal_risk': min(reference['values']),
                     'actual_max_moment_error': actual_error,
                     'coverage_observed': lower <= risk <= upper,
                     'moment_event_observed': actual_error <= learned['epsilon'],
                     'regret_bound_observed': risk - min(reference['values']) <= learned['conditional_regret_bound']})
    return serialize({'status': 'finite-sample calibration bridge; not a compounding result',
                      'physical_experiments_performed': False, 'cases': rows,
                      'learner_likelihood_table_access': False,
                      'independent_worlds': False,
                      'interpretation': 'two fixed source worlds, two calibration seeds each; standard empirical-moment estimator is an exact conventional tie'})


if __name__ == '__main__':
    result = run()
    out = ROOT / 'results/v4/calibration_results.json'
    out.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps([{'seed': r['seed'], 'risk': r['actual_risk'],
                       'interval': r['learned']['conditional_risk_interval'],
                       'moments_within_bound': r['moment_event_observed'],
                       'work': r['learned']['work']} for r in result['cases']], indent=2))
