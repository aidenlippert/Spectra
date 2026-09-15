"""Exploratory learned experiment-selection test with an exact strong baseline.

The learned object is a regression rule for two-shot action value. It learns
from charged solved design cases, not from supplied optimal-action code.
Physical models and their likelihoods remain supplied. This distinction is part
of the result, not a claim that these updates discover new physical laws.
"""
from dataclasses import asdict
from fractions import Fraction as F
from hashlib import sha256
import json
from math import log
from pathlib import Path
from time import perf_counter
import numpy as np
from experiments.v4_physics import make_family, likelihood_table, sample
from experiments.v4_planner import (branch, terminal_risk, eval_two_step_policy,
                                    verify_fixed_policy_certificate, new_work,
                                    compile_fixed_actions)
from experiments.v4_moments import signed_gram_plan
from experiments.v4_learner import CARTRegressor, NearestRecordRegressor

ROOT = Path(__file__).resolve().parents[1]
FEATURES = ('prior_sign_risk', 'one_shot_sign_risk', 'model_information_gain',
            'expected_max_model_mass', 'likelihood_variance', 'plus_probability',
            'outcome_imbalance', 'sign_information_gain', 'expected_posterior_purity',
            'channel_range', 'preparation_weight_fraction', 'prior_model_purity')


def entropy(probabilities):
    return -sum(float(p) * log(float(p)) for p in probabilities if p)


def add_work(a, b):
    out = dict(a)
    for key, value in b.items():
        out[key] = max(out.get(key, 0), value) if key == 'max_rational_bits' else out.get(key, 0) + value
    return out


def context(seed, n, prefix, denominator):
    rng = np.random.default_rng(seed)
    models, actions = make_family(n, seed, prefix)
    prep_v = F(int(rng.integers(45, 98)), denominator)
    readout_v = F(int(rng.integers(70, 99)), denominator)
    visibilities = tuple(prep_v ** sum(p != 'I' for p in a.prep) * readout_v for a in actions)
    table = likelihood_table(models, actions, visibilities)
    weights = [int(x) for x in rng.integers(1, 21, size=4)]
    prior = tuple(F(weights[i // 2], 2 * sum(weights)) for i in range(8))
    ticket = int(rng.integers(0, 2 * sum(weights)))
    truth = 0
    while ticket >= weights[truth // 2]:
        ticket -= weights[truth // 2]
        truth += 1
    history = []
    acquisition_work = new_work()
    for _ in range(2):
        a = int(rng.integers(0, min(2 * n, len(actions))))
        y = sample(models[truth], actions[a], 1, visibilities[a], rng)[0]
        probability, prior = branch(prior, table, a, y, acquisition_work)
        if not probability: raise AssertionError('source produced impossible observation')
        history.append({'action': a, 'outcome': y})
    public = {'models': [asdict(m) for m in models], 'actions': [asdict(a) for a in actions],
              'prior': [str(p) for p in prior], 'likelihoods': [[str(p) for p in row] for row in table],
              'visibilities': [str(v) for v in visibilities]}
    return {'seed': seed, 'n': n, 'prefix': prefix, 'models': models, 'actions': actions,
            'table': table, 'prior': prior, 'signs': tuple(m.target_sign for m in models),
            'visibilities': visibilities, 'history': history,
            'fingerprint': sha256(json.dumps(public, sort_keys=True).encode()).hexdigest(),
            'public': public, 'evaluator_hidden_model': truth,
            'acquisition_work': acquisition_work,
            'acquisition_cost': {'simulated_shots': 2, 'qubit_preparations': 2 * n,
                                 'source_gate_uses': 2 * len(models[0].gates)}}


def features(task):
    prior, table, signs = task['prior'], task['table'], task['signs']
    root_entropy = entropy(prior)
    sign_mass = [sum(p for p, s in zip(prior, signs) if s == z) for z in (-1, 1)]
    sign_entropy = entropy(sign_mass)
    work = new_work()
    rows = []
    for a, experiment in enumerate(task['actions']):
        risk, ent, sign_ent, max_mass, purity = 0., 0., 0., 0., 0.
        plus = 0.
        for y in (0, 1):
            probability, post = branch(prior, table, a, y, work)
            if y: plus = float(probability)
            if not probability: continue
            weight = float(probability)
            risk += weight * float(terminal_risk(post, signs))
            ent += weight * entropy(post)
            sign_ent += weight * entropy([sum(p for p, s in zip(post, signs) if s == z) for z in (-1, 1)])
            max_mass += weight * float(max(post))
            purity += weight * sum(float(p) ** 2 for p in post)
        variance = sum(float(p) * (float(table[i][a]) - plus) ** 2 for i, p in enumerate(prior))
        rows.append([float(terminal_risk(prior, signs)), risk, root_entropy - ent,
                     max_mass, variance, plus, abs(plus - .5), sign_entropy - sign_ent,
                     purity, float(max(row[a] for row in table) - min(row[a] for row in table)),
                     sum(p != 'I' for p in experiment.prep) / task['n'],
                     sum(float(p) ** 2 for p in prior)])
    return np.asarray(rows), work


def design_case(seed, n, prefix, denominator):
    task = context(seed, n, prefix, denominator)
    start = perf_counter()
    X, work = features(task)
    feature_seconds = perf_counter() - start
    start = perf_counter()
    teacher = signed_gram_plan(task['prior'], task['table'], task['signs'])
    teacher_seconds = perf_counter() - start
    return {'task': task, 'X': X, 'feature_work': work, 'feature_seconds': feature_seconds,
            'teacher': teacher, 'teacher_seconds': teacher_seconds,
            'y': np.asarray([float(v) for v in teacher['values']]) - X[:, 1]}


def collect(specs, denominator):
    return [design_case(seed, n, prefix, denominator) for seed, n, prefix in specs]


def fit_stages(first, second):
    X1, y1 = np.vstack([c['X'] for c in first]), np.concatenate([c['y'] for c in first])
    X2, y2 = np.vstack([c['X'] for c in second]), np.concatenate([c['y'] for c in second])
    all_x, all_y = np.vstack((X1, X2)), np.concatenate((y1, y2))
    models, fit_times = {}, {}
    def fit(name, X, y):
        start = perf_counter()
        models[name] = CARTRegressor(max_depth=4, min_leaf=8).fit(X, y)
        fit_times[name] = perf_counter() - start
    fit('first', X1, y1)
    fit('increment', X2, y2 - models['first'].predict(X2))
    fit('scratch', X2, y2)
    fit('scratch_increment', X2, y2 - models['scratch'].predict(X2))
    fit('batch', all_x, all_y)
    fit('batch_increment', all_x, all_y - models['batch'].predict(all_x))
    shuffled = np.random.default_rng(909).permutation(all_y)
    fit('shuffled', all_x, shuffled)
    start = perf_counter()
    models['nearest'] = NearestRecordRegressor().fit(all_x, all_y)
    fit_times['nearest'] = perf_counter() - start
    return models, fit_times


def proposals(case, models, exact_cache):
    X = case['X']; immediate = X[:, 1]
    goal = int(np.argmin(immediate)); info = int(np.argmax(X[:, 2]))
    purity = int(np.argmax(X[:, 8])); variance = int(np.argmax(X[:, 4]))
    data = {}
    def predict(name):
        start = perf_counter()
        value = models[name].predict(X)
        data[name] = {'seconds': perf_counter() - start, 'prediction_comparisons_or_distance_coordinates': models[name].predict_operations}
        return value
    first = predict('first'); increment = predict('increment')
    scratch = predict('scratch') + predict('scratch_increment')
    batch = predict('batch') + predict('batch_increment')
    nearest = predict('nearest'); shuffled = predict('shuffled')
    a1 = int(np.argmin(immediate + first))
    a2 = int(np.argmin(immediate + first + increment))
    # Erase action-dependent statistics, preserving only the global sign risk and
    # prior purity. This is an ablation, not a competitor granted equal features.
    erased = np.zeros_like(X); erased[:, 0] = X[:, 0]; erased[:, 11] = X[:, 11]
    start = perf_counter()
    erased_a = int(np.argmin(models['first'].predict(erased) + models['increment'].predict(erased)))
    erased_seconds = perf_counter() - start
    cached = exact_cache.get(case['task']['fingerprint'])
    rules = {
        'goal_greedy': [goal], 'entropy_greedy': [info],
        'fixed_portfolio': [goal, info, purity, variance],
        'raw_first': [a1], 'raw_cumulative': [a2],
        'frozen_learned': [goal, info, a1],
        'cumulative_learned': [goal, info, a1, a2],
        'scratch_second_stage': [goal, info, int(np.argmin(immediate + scratch)), variance],
        'batch_conventional_learner': [goal, info, int(np.argmin(immediate + batch)), variance],
        'nearest_record_retrieval': [goal, info, int(np.argmin(immediate + nearest)), variance],
        'exact_retrieval': [cached] if cached is not None else [goal, info, purity, variance],
        'shuffled_control': [goal, info, int(np.argmin(immediate + shuffled)), variance],
        'feature_erasure': [goal, info, erased_a, variance],
        'signed_gram_optimal': [case['teacher']['first_action']]}
    costs = {name: 0. for name in rules}
    for name in ('raw_first', 'frozen_learned'): costs[name] = data['first']['seconds']
    for name in ('raw_cumulative', 'cumulative_learned'): costs[name] = data['first']['seconds'] + data['increment']['seconds']
    costs['scratch_second_stage'] = data['scratch']['seconds'] + data['scratch_increment']['seconds']
    costs['batch_conventional_learner'] = data['batch']['seconds'] + data['batch_increment']['seconds']
    costs['nearest_record_retrieval'] = data['nearest']['seconds']
    costs['shuffled_control'] = data['shuffled']['seconds']
    costs['feature_erasure'] = erased_seconds
    return {k: sorted(set(v)) for k, v in rules.items()}, costs, data, cached is not None


def evaluate_case(case, models, cache):
    task = case['task']; p, L, signs = task['prior'], task['table'], task['signs']
    choices, inference_seconds, inference_detail, cache_hit = proposals(case, models, cache)
    outputs = {}
    for name, candidates in choices.items():
        start = perf_counter()
        generation, verification = new_work(), new_work()
        certificates = []
        for action in candidates:
            if name == 'signed_gram_optimal':
                cert = compile_fixed_actions(p, L, signs, action,
                         case['teacher']['child_actions'][action], generation)
            else:
                cert = eval_two_step_policy(p, L, signs, action, counter=generation)
            if not verify_fixed_policy_certificate(cert, p, L, signs, verification):
                raise AssertionError('independent policy risk certificate rejected')
            if F(cert['risk']) != case['teacher']['values'][action]:
                raise AssertionError('signed moment and posterior policy values disagree')
            certificates.append(cert)
        selected = min(certificates, key=lambda c: (F(c['risk']), c['first_action']))
        measured_seconds = perf_counter() - start
        cost = {'generation': generation, 'independent_verification': verification,
                'candidate_count': len(candidates), 'measured_policy_seconds': measured_seconds,
                'measured_inference_seconds': inference_seconds[name]}
        if name == 'signed_gram_optimal':
            cost['teacher'] = case['teacher']['work']
            cost['measured_total_seconds'] = measured_seconds + case['teacher_seconds']
        else:
            cost['feature_extraction'] = case['feature_work']
            cost['measured_total_seconds'] = measured_seconds + inference_seconds[name] + case['feature_seconds']
        outputs[name] = {'risk': selected['risk'], 'regret_vs_optimal': str(F(selected['risk']) - min(case['teacher']['values'])),
                         'certificate': selected, 'all_candidate_certificates': certificates,
                         'work': cost, 'target_shots': 2, 'target_qubit_preparations': 2 * task['n'],
                         'target_source_gate_uses': 2 * len(task['models'][0].gates)}
    return {'seed': task['seed'], 'n': task['n'], 'prefix': task['prefix'], 'fingerprint': task['fingerprint'],
            'context': task['public'], 'calibration_history': task['history'],
            'acquisition_cost': task['acquisition_cost'], 'exact_retrieval_hit': cache_hit,
            'methods': outputs, 'inference_detail': inference_detail}


def serial_case(case):
    task = case['task']
    return {'seed': task['seed'], 'fingerprint': task['fingerprint'], 'n': task['n'], 'prefix': task['prefix'],
            'context': task['public'], 'calibration_history': task['history'],
            'acquisition_cost': task['acquisition_cost'],
            'features': case['X'].tolist(), 'teacher_residual_targets': case['y'].tolist(),
            'teacher_action_values': [str(v) for v in case['teacher']['values']],
            'teacher_first_action': case['teacher']['first_action'], 'teacher_work': case['teacher']['work'],
            'teacher_seconds': case['teacher_seconds'], 'feature_work': case['feature_work']}


def run(small=False):
    started = perf_counter()
    count1, count2, count_test = (6, 8, 6) if small else (24, 32, 32)
    specs1 = [(101 + i, 2 + i % 2, False) for i in range(count1)]
    specs2 = [(401 + i, 2 + i % 2, bool(i % 2)) for i in range(count2)]
    specs_test = [(1001 + i, 3 + i % 2, bool((i // 2) % 2)) for i in range(count_test)]
    first, second = collect(specs1, 101), collect(specs2, 101)
    models, fit_times = fit_stages(first, second)
    cache = {c['task']['fingerprint']: c['teacher']['first_action'] for c in first + second}
    test = collect(specs_test, 103)
    rows = [evaluate_case(case, models, cache) for case in test]
    summary = {}
    for method in rows[0]['methods']:
        risks = [F(row['methods'][method]['risk']) for row in rows]
        regrets = [F(row['methods'][method]['regret_vs_optimal']) for row in rows]
        summary[method] = {'mean_exact_risk': float(sum(risks) / len(risks)),
                           'mean_exact_regret': float(sum(regrets) / len(regrets)),
                           'worst_exact_regret': float(max(regrets)),
                           'exact_optimum_matches': sum(r == 0 for r in regrets),
                           'mean_measured_cpu_seconds': sum(row['methods'][method]['work']['measured_total_seconds'] for row in rows) / len(rows),
                           'mean_candidate_count': sum(row['methods'][method]['work']['candidate_count'] for row in rows) / len(rows)}
    combined = first + second
    return {'status': 'exploratory_protocol; not an externally blinded preregistration',
            'protocol': {'first_stage_tasks': specs1, 'second_stage_tasks': specs2,
                         'test_tasks': specs_test, 'training_noise_denominator': 101,
                         'test_noise_denominator': 103, 'feature_names': FEATURES,
                         'model_class': 'one/two coupled Pauli Clifford gates; supplied candidate set and calibrated noise'},
            'training_cases': [serial_case(c) for c in combined],
            'models': {k: v.to_dict() for k, v in models.items() if k != 'nearest'},
            'fit_work': {k: {'fit_seconds': fit_times[k],
                              'split_candidates': getattr(v, 'split_candidates', None),
                              'sse_evaluations': getattr(v, 'sse_evaluations', None),
                              'featurevalues_scanned': getattr(v, 'featurevalues_scanned', None)} for k, v in models.items()},
            'training_cost': {'simulated_calibration_shots': sum(c['task']['acquisition_cost']['simulated_shots'] for c in combined),
                              'simulated_qubit_preparations': sum(c['task']['acquisition_cost']['qubit_preparations'] for c in combined),
                              'teacher_wL_products': sum(c['teacher']['work']['wL_products'] for c in combined),
                              'teacher_gram_products': sum(c['teacher']['work']['G_products'] for c in combined),
                              'teacher_cpu_seconds': sum(c['teacher_seconds'] for c in combined)},
            'evaluation': rows, 'summary': summary, 'exact_cache_hits': sum(r['exact_retrieval_hit'] for r in rows),
            'runtime_seconds': perf_counter() - started, 'physical_experiments_performed': False,
            'interpretation': 'data-learned planning-value rules; physics and Bayesian likelihood model supplied; measured gates must be evaluated before any compounding claim'}


if __name__ == '__main__':
    result = run()
    out = ROOT / 'results/v4'; out.mkdir(parents=True, exist_ok=True)
    (out / 'policy_results.json').write_text(json.dumps(result, indent=2) + '\n')
    (out / 'summary.json').write_text(json.dumps({k: result[k] for k in ('status', 'protocol', 'training_cost', 'fit_work', 'summary', 'exact_cache_hits', 'runtime_seconds', 'interpretation')}, indent=2) + '\n')
    print(json.dumps(result['summary'], indent=2))
