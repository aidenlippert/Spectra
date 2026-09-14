"""Verify persisted V4 source tables, policy risks and calibration guarantees."""
from pathlib import Path
from fractions import Fraction as F
from hashlib import sha256
import json
import re
from experiments.v4_planner import verify_fixed_policy_certificate
from experiments.v4_moments import signed_gram_plan
from experiments.v4_physics import Model, Experiment, likelihood_table
from experiments.v4_calibration import test_world, exact_policy_risk, sample_budget

ROOT = Path(__file__).resolve().parents[1]


def verify():
    r = json.loads((ROOT / 'results/v4/policy_results.json').read_text())
    checked = selected = 0
    for row in r['evaluation']:
        c = row['context']
        p = tuple(map(F, c['prior']))
        L = tuple(tuple(map(F, x)) for x in c['likelihoods'])
        models = tuple(Model(tuple((str(g[0]), int(g[1])) for g in item['gates']), item['target_sign']) for item in c['models'])
        actions = tuple(Experiment(**item) for item in c['actions'])
        assert likelihood_table(models, actions, tuple(map(F, c['visibilities']))) == L
        signs = tuple(m.target_sign for m in models)
        optimum = min(signed_gram_plan(p, L, signs)['values'])
        for item in row['methods'].values():
            assert verify_fixed_policy_certificate(item['certificate'], p, L, signs)
            assert F(item['regret_vs_optimal']) == F(item['risk']) - optimum
            for cert in item['all_candidate_certificates']:
                assert verify_fixed_policy_certificate(cert, p, L, signs)
                checked += 1
            selected += 1
            assert item['target_qubit_preparations'] == 2 * row['n']
        assert row['methods']['cumulative_learned']['risk'] == row['methods']['frozen_learned']['risk']
    train = {c['fingerprint'] for c in r['training_cases']}
    assert all(row['fingerprint'] not in train for row in r['evaluation'])
    cal = json.loads((ROOT / 'results/v4/calibration_results.json').read_text())
    for row in cal['cases']:
        p, L, signs = test_world(row['right_world'])
        policy = row['learned']
        actual = exact_policy_risk(policy, p, L, signs)
        assert actual == F(row['actual_risk']) == F(9, 50)
        assert actual == min(signed_gram_plan(p, L, signs)['values'])
        assert F(policy['conditional_risk_interval'][0]) <= actual <= F(policy['conditional_risk_interval'][1])
        count, n, power = sample_budget(4, F(policy['epsilon']), F(policy['delta']))
        assert policy['work']['labelled_episodes'] == count * n
        assert policy['work']['visible_readouts'] == 24 * n
        assert F(policy['union_failure_upper_bound']) == 2 * count * F(1, 2 ** power) <= F(policy['delta'])
    log = (ROOT / 'results/v4/test_log.txt').read_text()
    tests = int(re.search(r'Ran (\d+) tests', log).group(1))
    assert log.rstrip().endswith('OK')
    paths = sorted(list(ROOT.glob('experiments/v4_*.py')) + list(ROOT.glob('tests/test_v4_*.py')))
    receipt = {'date': '2026-09-10', 'test_count': tests, 'tests_passed': True,
               'saved_selected_certificates_verified': selected,
               'saved_candidate_certificates_verified': checked,
               'source_likelihood_tables_recomputed': len(r['evaluation']),
               'calibration_policy_risks_and_counts_verified': len(cal['cases']),
               'train_test_exact_context_overlap': 0,
               'second_learned_increment_positive': False,
               'calibration_learner_likelihood_table_access': False,
               'physical_experiments_performed': False, 'larger_goal_complete': False,
               'hashes': {str(p.relative_to(ROOT)): sha256(p.read_bytes()).hexdigest() for p in paths},
               'result_hashes': {name: sha256((ROOT / 'results/v4' / name).read_bytes()).hexdigest()
                                 for name in ('policy_results.json', 'calibration_results.json')},
               'limits': ['V4 policy distillation has supplied physics and calibrated noise.',
                          'The calibration bridge uses supplied target labels and moment vocabulary.',
                          'Source conditions and physical preparation costs remain model assumptions.',
                          'Operation counters and timings are not complete bit-complexity bounds.',
                          'No population generalization or broad scientific-compounding claim.']}
    (ROOT / 'results/v4/verification_receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    return receipt


if __name__ == '__main__':
    result = verify()
    print(json.dumps({k: v for k, v in result.items() if k not in ('hashes', 'limits', 'result_hashes')}, indent=2))
