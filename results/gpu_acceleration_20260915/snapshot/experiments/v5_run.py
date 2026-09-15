"""Two sequential, costed transfers in unknown dense coupled Gaussian matter.

Verification contexts are simulated. No selected eigenmode, target sign or
random seed is passed to a learner. All sources have the same charged reset
and finite equilibration period; exact-cache task keys never repeat.
"""
from dataclasses import asdict
from fractions import Fraction as F
from pathlib import Path
from time import perf_counter
import json
import numpy as np
from experiments.v5_core import (Contract, LearnedPrefix, acquire_next, make_source,
    replay_records, resource_lower_bounds, PREPARATION_TIME)
from experiments.v5_frame import null_probe
from experiments.v5_bayes import (prior_empty, extend_belief, update_belief,
    choose_probe, current_target_decision)

ROOT = Path(__file__).resolve().parents[1]
METHODS = ('cumulative', 'frozen_after_first', 'scratch', 'exact_lookup',
           'structured_retrieval', 'stateful_bayes', 'fixed_one_read')


def serialize(value):
    if isinstance(value, F): return str(value)
    if isinstance(value, dict): return {str(k): serialize(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)): return [serialize(x) for x in value]
    return value


def world(seed, branches=4, targets=4):
    rng = np.random.default_rng(seed)
    a = int(rng.choice((-1, 1)))
    alpha_numerators = rng.choice(np.arange(1001, 2000), size=branches * (targets + 1), replace=False)
    cursor = 0
    parent = f'component-{int(rng.integers(10**8, 10**9))}'
    out = {'setup': {'signs': (a,), 'contract': Contract((parent,), F(1)),
                     'source_seed': int(rng.integers(1, 10**8))}, 'branches': []}
    for j in range(branches):
        b = int(rng.choice((-1, 1)))
        child = f'component-{int(rng.integers(10**8, 10**9))}'
        lineage = (parent, child)
        base = {'signs': (a, b), 'contract': Contract(lineage, F(int(alpha_numerators[cursor]), 1000)),
                'source_seed': int(rng.integers(1, 10**8)), 'targets': []}
        cursor += 1
        for k in range(targets):
            c = int(rng.choice((-1, 1)))
            grandchild = f'component-{int(rng.integers(10**8, 10**9))}'
            base['targets'].append({'signs': (a, b, c),
                'contract': Contract(lineage + (grandchild,), F(int(alpha_numerators[cursor]), 1000)),
                'source_seed': int(rng.integers(1, 10**8))})
            cursor += 1
        out['branches'].append(base)
    return out


def run_world(data, method):
    if method not in METHODS: raise ValueError('unknown comparison method')
    started = perf_counter()
    output = {'method': method, 'tasks': [], 'source_ledgers': [],
              'source_records': [],
              'likelihood_evaluations': 0, 'record_replay_updates': 0,
              'exact_cache_hits': 0, 'calibration_label_reads': 0}
    cache = {}
    def source(item):
        read, ledger, records = make_source(item['signs'], item['contract'], item['source_seed'])
        output['source_ledgers'].append(ledger)
        output['source_records'].append({'lineage': item['contract'].lineage,
                                         'records': records})
        return read
    def note(item, estimate, learned=None):
        # Truth is confined to evaluation, after the selected policy has run.
        record = {'lineage': item['contract'].lineage, 'alpha': str(item['contract'].alpha),
                  'stage': len(item['contract'].lineage), 'estimate': estimate,
                  'evaluator_hidden_prefix': item['signs'],
                  'evaluator_target': item['signs'][-1], 'correct': estimate == item['signs'][-1]}
        if learned is not None:
            record['learned'] = asdict(learned)
        output['tasks'].append(record)
    def from_scratch(item):
        read = source(item); learned = LearnedPrefix()
        for level in range(1, len(item['contract'].lineage) + 1):
            contract = Contract(item['contract'].lineage[:level], item['contract'].alpha)
            learned = acquire_next(read, contract, learned)
        return learned
    def bayes_read(item, belief):
        read = source(item); probe = choose_probe(belief)
        value = read(len(item['contract'].lineage), probe)
        updated = update_belief(extend_belief(belief), len(item['contract'].lineage), item['contract'].alpha, probe, value)
        output['likelihood_evaluations'] += len(updated)
        return updated

    initial = None; initial_belief = None
    if method in ('cumulative', 'frozen_after_first', 'structured_retrieval'):
        initial = acquire_next(source(data['setup']), data['setup']['contract'])
        output['initial_artifact'] = asdict(initial)
    elif method == 'stateful_bayes':
        initial_belief = bayes_read(data['setup'], prior_empty())
    for branch in data['branches']:
        learned_b = belief_b = None
        if method in ('cumulative', 'frozen_after_first', 'structured_retrieval'):
            previous = replay_records(initial.records) if method == 'structured_retrieval' else initial
            if method == 'structured_retrieval': output['record_replay_updates'] += len(initial.records)
            learned_b = acquire_next(source(branch), branch['contract'], previous)
            note(branch, learned_b.signs[-1], learned_b)
        elif method == 'stateful_bayes':
            belief_b = bayes_read(branch, initial_belief)
            note(branch, current_target_decision(belief_b))
        elif method in ('scratch', 'exact_lookup'):
            learned_b = from_scratch(branch)
            note(branch, learned_b.signs[-1], learned_b)
            cache[(branch['contract'].lineage, branch['contract'].alpha)] = learned_b.signs[-1]
        else:
            y = source(branch)(2, null_probe((1,)))
            note(branch, 1 - 2 * y)
        for target in branch['targets']:
            if method in ('cumulative', 'structured_retrieval'):
                previous = replay_records(learned_b.records) if method == 'structured_retrieval' else learned_b
                if method == 'structured_retrieval': output['record_replay_updates'] += len(learned_b.records)
                learned = acquire_next(source(target), target['contract'], previous)
                note(target, learned.signs[-1], learned)
            elif method == 'frozen_after_first':
                read = source(target)
                transient_b = acquire_next(read, Contract(target['contract'].lineage[:2], target['contract'].alpha), initial)
                learned = acquire_next(read, target['contract'], transient_b)
                note(target, learned.signs[-1], learned)
            elif method == 'stateful_bayes':
                posterior = bayes_read(target, belief_b)
                note(target, current_target_decision(posterior))
                # Retain posterior information about the shared prefix after
                # integrating out this target's fresh third-mode sign.
                from experiments.v5_bayes import Belief
                belief_b = Belief({h: sum(w for full, w in posterior.items() if full[:-1] == h)
                                   for h in belief_b}, posterior.workcount)
            elif method in ('scratch', 'exact_lookup'):
                key = (target['contract'].lineage, target['contract'].alpha)
                if method == 'exact_lookup' and key in cache:
                    output['exact_cache_hits'] += 1; note(target, cache[key])
                else:
                    learned = from_scratch(target)
                    note(target, learned.signs[-1], learned)
                    cache[key] = learned.signs[-1]
            else:
                y = source(target)(3, null_probe((1, 1)))
                note(target, 1 - 2 * y)
    output['total_scalar_reads'] = sum(x['scalar_readouts'] for x in output['source_ledgers'])
    output['total_source_preparations'] = sum(x['source_preparations'] for x in output['source_ledgers'])
    output['total_relaxation_time_units'] = sum(x['relaxation_time_units'] for x in output['source_ledgers'])
    output['max_record_coordinate_bits'] = max((max(F(x).numerator.bit_length(), F(x).denominator.bit_length())
                                              for row in output['tasks'] if 'learned' in row
                                              for record in row['learned']['records'] for x in record['probe']), default=0)
    output['cpu_seconds'] = perf_counter() - started
    return output


def run(count=64, branches=4, targets=4):
    if not 1 <= count <= 128 or branches != 4 or targets != 4:
        raise ValueError('bounded verification design requires four-by-four transfer')
    results = []
    for index in range(count):
        seed = 41001 + index
        data = world(seed, branches, targets)
        results.append({'evaluator_seed': seed, 'methods': {method: run_world(data, method) for method in METHODS}})
    summary = {}
    for method in METHODS:
        rows = [r['methods'][method] for r in results]
        summary[method] = {'worlds': count, 'source_reads': sum(x['total_scalar_reads'] for x in rows),
                          'mean_source_reads_per_world': sum(x['total_scalar_reads'] for x in rows)/count,
                          'source_preparations': sum(x['total_source_preparations'] for x in rows),
                          'relaxation_time_units': sum(x['total_relaxation_time_units'] for x in rows),
                          'exact_cache_hits': sum(x['exact_cache_hits'] for x in rows),
                          'likelihood_evaluations': sum(x['likelihood_evaluations'] for x in rows),
                          'record_replay_updates': sum(x['record_replay_updates'] for x in rows),
                          'mean_cpu_seconds_per_world': sum(x['cpu_seconds'] for x in rows)/count}
        for stage in (2, 3):
            tasks = [item for row in rows for item in row['tasks'] if item['stage'] == stage]
            summary[method][f'stage_{stage}_tasks'] = len(tasks)
            summary[method][f'stage_{stage}_errors'] = sum(not x['correct'] for x in tasks)
            summary[method][f'stage_{stage}_error_fraction'] = sum(not x['correct'] for x in tasks)/len(tasks)
    return serialize({'scope': 'bounded coupled-Gaussian representation acquisition and two costed transfers',
                      'physical_experiments_performed': False,
                      'protocol': {'worlds': count, 'branches_per_world': branches,
                                   'third_stage_targets_per_branch': targets,
                                   'evaluator_seeds': [41001 + i for i in range(count)],
                                   'coupling_range': '[1,2] times (1e14,1e10,1e6)',
                                   'control_error_norm': '1/1000000000', 'readout_error': '1/1000',
                                   'source_preparation_time': PREPARATION_TIME,
                                   'same_normal_random_numbers_across_methods': True,
                                   'learner_access': 'public contract and binary scalar readouts; no source labels'},
                      'analytical_bounds': resource_lower_bounds(),
                      'summary': summary, 'world_results': results})


if __name__ == '__main__':
    result = run()
    directory = ROOT / 'results/v5'; directory.mkdir(parents=True, exist_ok=True)
    (directory / 'transfer_results.json').write_text(json.dumps(result, indent=2) + '\n')
    (directory / 'summary.json').write_text(json.dumps({k: result[k] for k in ('scope', 'protocol', 'analytical_bounds', 'summary')}, indent=2) + '\n')
    print(json.dumps(result['summary'], indent=2))
