"""Auditable, finite demonstration of sequential response-law acquisition.

The learner receives calibration records. Hidden masks are used only by the
source and evaluator. Every scientific guarantee is conditional on the supplied
GF(2) law family and declared quantum source/measurement menu.
"""
from __future__ import annotations
from dataclasses import asdict
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
import json
from math import log, sqrt
from pathlib import Path
from time import perf_counter
import numpy as np
from experiments.v2_laws import CalibrationRecord, learn, predict, required_shots, verify_certificate
from experiments.v2_physics import World, born_probability
from experiments.v2_policy import certificate, verify

ROOT = Path(__file__).resolve().parents[1]
METHODS = ('frozen', 'A_only', 'B_only', 'AB', 'exact_lookup',
           'structured_retrieval', 'conventional_GF2')


def contexts(d):
    training = tuple(tuple(int(i == j) for i in range(d)) for j in range(d))
    heldout = tuple(x for x in product((0, 1), repeat=d) if any(x) and x not in training)
    return training, heldout


def compiled_b_error(visibility):
    """Prove both transition probabilities of the randomized compiled channel."""
    if not isinstance(visibility, F) or not 0 < visibility <= 1:
        raise ValueError('positive exact visibility required')
    r = visibility / (2 + visibility)
    p1_if_b0 = (1 - visibility) / 2 + (1 + visibility) / 2 * r
    p1_if_b1 = (1 + r) / 2
    eta = 1 / (2 + visibility)
    if p1_if_b0 != eta or p1_if_b1 != 1 - eta:
        raise AssertionError('compiled channel identity failed')
    return eta, r


def conventional_fit(records):
    """Independent exhaustive linear-model baseline using identical raw records.

    Chooses the unique mask agreeing with all majority labels. For these square
    independent designs it is also the joint maximum-likelihood mask under BSC
    noise. Its search cost is explicitly counted; this is not an optimal runtime
    baseline. Gaussian elimination can implement the same inference as cheaply
    as the principal learner.
    """
    d = len(records[0].context)
    labels = tuple(int(sum(r.outcomes) > len(r.outcomes) // 2) for r in records)
    candidates = []
    bit_products = 0
    for mask in product((0, 1), repeat=d):
        agrees = True
        for rec, label in zip(records, labels):
            bit_products += d
            if sum(a * b for a, b in zip(mask, rec.context)) % 2 != label:
                agrees = False
                break
        if agrees:
            candidates.append(mask)
    if len(candidates) != 1:
        raise ValueError('conventional fit does not identify one mask')
    return candidates[0], bit_products


def acquire(a_sensor, reference_sensor, d, visibility, delta_a, delta_b, rng):
    """The only acquisition interface is sampled-data callbacks, never World."""
    training, heldout = contexts(d)
    if not heldout:
        raise ValueError('need unseen contexts')
    eta_a = F(1, 10)
    eta_b, postprocess = compiled_b_error(visibility)
    ka = required_shots(d, eta_a, delta_a)
    kb = required_shots(d, eta_b, delta_b)
    a_records = tuple(CalibrationRecord(x, tuple(a_sensor(x, ka))) for x in training)
    a_result = learn(a_records, eta_a, delta_a)
    if not verify_certificate(a_records, a_result, eta_a, delta_a):
        raise ValueError('A acquisition lacks a valid model-conditional certificate')
    reference_x = heldout[0]
    compiled_axis = predict(a_result.mask, reference_x)
    b_records, raw_reference, coins = [], [], []
    for z in training:
        outcomes = tuple(reference_sensor(reference_x, z, (compiled_axis, 0), kb))
        # A plus result becomes bit 1 with probability r; a minus is always 1.
        coin = tuple(int(t < postprocess.numerator) for t in
                     rng.integers(0, postprocess.denominator, size=kb))
        labels = tuple(1 if y == -1 else t for y, t in zip(outcomes, coin))
        b_records.append(CalibrationRecord(z, labels))
        raw_reference.append(outcomes)
        coins.append(coin)
    b_records = tuple(b_records)
    b_result = learn(b_records, eta_b, delta_b)
    if not verify_certificate(b_records, b_result, eta_b, delta_b):
        raise ValueError('B acquisition lacks a valid conditional certificate')
    a_conventional, a_ops = conventional_fit(a_records)
    b_conventional, b_ops = conventional_fit(b_records)
    return {
        'a': a_result, 'b': b_result,
        'a_records': a_records, 'b_records': b_records,
        'raw_reference': raw_reference, 'postprocessing_coins': coins,
        'reference_x': reference_x, 'compiled_axis': compiled_axis,
        'eta_b': eta_b, 'postprocess': postprocess,
        'conventional': (a_conventional, b_conventional),
        'cost': {
            'a_single_qubit_preparations_and_shots': d * ka,
            'b_two_qubit_preparations_and_shots': d * kb,
            'total_acquisition_shots': d * (ka + kb),
            'acquisition_qubits_prepared': d * ka + 2 * d * kb,
            'postprocessing_coin_draws': d * kb,
            'learned_rule_storage_bits': 2 * d,
            'retained_record_outcome_bits': d * (ka + 3 * kb),
            'context_bits_in_calibration_records': 2 * d * d,
            'learner_gaussian_elimination_xors': a_result.xor_operations + b_result.xor_operations,
            'verifier_gaussian_elimination_xors': a_result.xor_operations + b_result.xor_operations,
            'majority_bits_scanned_per_learner_or_verifier_pass': d * (ka + kb),
            'conventional_enumeration_bit_products': a_ops + b_ops,
            'shots_per_A_context': ka, 'shots_per_B_context': kb,
            'certificate_arithmetic': 'exact binomial-tail integers; wall time reported separately',
        },
    }


def exact_expected_risk(method, x, z, truth_a, truth_b, acquired):
    """Evaluator only: integrate random guesses and measurement exactly."""
    ar, br = acquired['a'].mask, acquired['b'].mask
    ai = (predict(ar, x),) if method in ('A_only', 'AB', 'structured_retrieval') else (0, 1)
    bj = (predict(br, z),) if method in ('B_only', 'AB', 'structured_retrieval') else (0, 1)
    if method == 'conventional_GF2':
        ai, bj = (predict(acquired['conventional'][0], x),), (predict(acquired['conventional'][1], z),)
    # On held-out contexts an exact input-answer table has no entry.
    total = sum(F(1, 2) - F(1, 4) * (i == truth_a) * (j == truth_b)
                for i in ai for j in bj)
    return total / (len(ai) * len(bj))


def evaluate(world, acquired, heldout, trials, rng, family_failure, comparisons):
    counts = {m: 0 for m in METHODS}
    exact = {m: F(0) for m in METHODS}
    conventional_ties = True
    lookup = ({r.context: int(sum(r.outcomes) > len(r.outcomes) // 2)
               for r in acquired['a_records']},
              {r.context: int(sum(r.outcomes) > len(r.outcomes) // 2)
               for r in acquired['b_records']})
    lookup_hits = 0
    ar, br = acquired['a'].mask, acquired['b'].mask
    # Vectorized independent task sampling; a common outcome random number makes
    # within-task contrasts paired. Tasks remain independent conditional on laws.
    indices = rng.integers(0, len(heldout), size=(trials, 2))
    guesses = rng.integers(0, 2, size=(trials, 2))
    signs = rng.integers(0, 2, size=trials) * 2 - 1
    outcome_draws = rng.integers(0, 4, size=trials)
    pairs = set()
    for t in range(trials):
        x, z = heldout[int(indices[t, 0])], heldout[int(indices[t, 1])]
        pairs.add((x, z))
        lookup_hits += int(x in lookup[0]) + int(z in lookup[1])
        # Evaluator truth is never passed to acquisition or action selection.
        ta, tb = predict(world.mask_a, x), predict(world.mask_b, z)
        a, b = predict(ar, x), predict(br, z)
        ca, cb = predict(acquired['conventional'][0], x), predict(acquired['conventional'][1], z)
        conventional_ties &= (a, b) == (ca, cb)
        ga, gb = (int(v) for v in guesses[t])
        actions = {'frozen': (ga, gb), 'A_only': (a, gb), 'B_only': (ga, b),
                   'AB': (a, b), 'exact_lookup': (lookup[0].get(x, ga), lookup[1].get(z, gb)),
                   'structured_retrieval': (a, b), 'conventional_GF2': (ca, cb)}
        for method, action in actions.items():
            p = born_probability(ta, tb, int(signs[t]), action)
            outcome = 1 if F(int(outcome_draws[t]), 4) < p else -1
            counts[method] += int(outcome != signs[t])
    # Exhaustive unseen contexts make the conditional expected-risk receipt exact.
    for x, z in product(heldout, repeat=2):
        ta, tb = predict(world.mask_a, x), predict(world.mask_b, z)
        for method in METHODS:
            exact[method] += exact_expected_risk(method, x, z, ta, tb, acquired)
    exact = {m: v / len(heldout) ** 2 for m, v in exact.items()}
    radius = sqrt(log(2 * comparisons / float(family_failure)) / (2 * trials))
    rows = {m: {'errors': counts[m], 'trials': trials,
                'empirical_error': counts[m] / trials,
                'simultaneous_hoeffding_interval': [max(0, counts[m] / trials - radius),
                                                   min(1, counts[m] / trials + radius)],
                'exact_conditional_risk': str(exact[m])} for m in METHODS}
    return {'methods': rows, 'exact_complementarity': str(exact['A_only'] + exact['B_only'] - exact['frozen'] - exact['AB']),
            'unique_target_pairs_sampled': len(pairs),
            'unseen_target_pairs_exhaustively_checked': len(heldout) ** 2,
            'exact_lookup_hits': lookup_hits, 'conventional_action_tie': conventional_ties,
            'shared_target_shots_per_method': trials,
            'total_simulated_method_shots': trials * len(METHODS),
            'paired_randomness': True,
            'confidence_scope': 'conditional on fixed acquired laws; union over all worlds and methods'}


def serial_acquisition(acq):
    out = {k: v for k, v in acq.items() if k not in ('a', 'b', 'a_records', 'b_records', 'eta_b', 'postprocess')}
    out.update(a=asdict(acq['a']), b=asdict(acq['b']),
               a_records=[asdict(r) for r in acq['a_records']],
               b_records=[asdict(r) for r in acq['b_records']],
               eta_b=str(acq['eta_b']), postprocess=str(acq['postprocess']))
    return out


def run(protocol=None):
    started = perf_counter()
    if protocol is None:
        protocol = json.loads((ROOT / 'research/v2/protocol.json').read_text())
    # The theorem and evaluator here are for this fixed visibility and dimensions.
    if protocol['visibility'] != '1/2' or protocol['a_channel_error'] != '1/10':
        raise ValueError('unsupported protocol channel')
    d, trials = protocol['dimension'], protocol['evaluation_trials_per_world']
    if type(d) is not int or not 2 <= d <= 8 or type(trials) is not int or not 1 <= trials <= 100000:
        raise ValueError('protocol outside finite resource budget')
    seeds = protocol['world_seeds']
    if not 1 <= len(seeds) <= 16 or any(type(s) is not int or s < 0 for s in seeds):
        raise ValueError('invalid seed budget')
    da, db = F(protocol['delta_a']), F(protocol['delta_b_conditional_on_a'])
    pc_start = perf_counter()
    certs = {f'{known}_h{h}': certificate(known, h)
             for known, h in [('none', 1), ('A', 1), ('B', 1), ('AB', 1), ('none', 4)]}
    if not all(verify(c) for c in certs.values()):
        raise AssertionError('optimal-policy certificate failed')
    policy_seconds = perf_counter() - pc_start
    ideal_j = F(certs['A_h1']['value']) + F(certs['B_h1']['value']) - F(certs['none_h1']['value']) - F(certs['AB_h1']['value'])
    if ideal_j != F(1, 16) or F(certs['none_h4']['value']) <= F(3, 10):
        raise AssertionError('theorem obligation failed')
    rows = []
    for seed in seeds:
        source_seed, noise_seed, post_seed, eval_seed = np.random.SeedSequence(seed).spawn(4)
        source_rng, noise_rng, post_rng, eval_rng = (np.random.default_rng(s) for s in
                                                   (source_seed, noise_seed, post_seed, eval_seed))
        ma, mb = (tuple(int(b) for b in source_rng.integers(0, 2, size=d)) for _ in range(2))
        world = World(ma, mb)
        start = perf_counter()
        acq = acquire(lambda x, n: world.calibration('a', x, n, noise_rng),
                      lambda x, z, action, n: world.target(x, z, 1, action, n, noise_rng),
                      d, F(1, 2), da, db, post_rng)
        acquisition_seconds = perf_counter() - start
        _, heldout = contexts(d)
        start = perf_counter()
        ev = evaluate(world, acq, heldout, trials, eval_rng,
                      F(protocol['evaluation_family_failure']), len(METHODS) * len(seeds))
        eval_seconds = perf_counter() - start
        total_delta = F(acq['a'].exact_error_bound) + F(acq['b'].exact_error_bound)
        acquisition_shots = acq['cost']['total_acquisition_shots']
        rows.append({'seed': seed, 'evaluator_truth': {'mask_a': ma, 'mask_b': mb},
                     'acquisition': serial_acquisition(acq), 'evaluation': ev,
                     'both_masks_recovered': acq['a'].mask == ma and acq['b'].mask == mb,
                     'guarantees': {'acquisition_failure_union_bound': str(total_delta),
                                    'AB_unconditional_error_upper_bound': str(F(1, 4) + total_delta / 4),
                                    'complementarity_unconditional_lower_bound': str(F(1, 16) - total_delta / 4)},
                     'fixed_horizon_budget': {'target_error_threshold': '3/10',
                                               'reset_target_only_minimum_shots_per_task_lower_bound': 5,
                                               'learner_acquisition_plus_target_shots': acquisition_shots + trials,
                                               'reset_target_only_shot_budget_lower_bound': 5 * trials,
                                               'strict_shot_break_even_tasks': acquisition_shots // 4 + 1,
                                               'scope': 'committed fixed-horizon budgets; not expected stopping times or globally learning baseline'},
                     'runtime_seconds': {'acquisition_with_independent_verification_and_baseline': acquisition_seconds,
                                         'evaluation_and_exhaustive_risk_check': eval_seconds}})
    return {'protocol': protocol,
            'protocol_sha256': sha256(json.dumps(protocol, sort_keys=True).encode()).hexdigest(),
            'theorem': {'one_shot_risks': {k: certs[f'{k}_h1']['value'] for k in ('none', 'A', 'B', 'AB')},
                        'ideal_complementarity': str(ideal_j),
                        'optimal_reset_four_shot_risk': certs['none_h4']['value'],
                        'certificate_nodes_checked': sum(len(c['nodes']) for c in certs.values())},
            'policy_certificates': certs, 'worlds': rows,
            'runtime_seconds': {'policy_generation_and_checking': policy_seconds,
                                'total': perf_counter() - started},
            'claim': 'finite model-conditional sequential acquisition and complementary transfer; conventional structured inference ties',
            'physical_experiments_performed': False}


def main():
    result = run()
    out = ROOT / 'results/v2'
    out.mkdir(parents=True, exist_ok=True)
    (out / 'compounding_results.json').write_text(json.dumps(result, indent=2) + '\n')
    summary = {k: v for k, v in result.items() if k not in ('policy_certificates', 'worlds')}
    summary['worlds'] = [{k: v for k, v in w.items() if k != 'acquisition'} | {'acquisition_cost': w['acquisition']['cost']}
                         for w in result['worlds']]
    (out / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
