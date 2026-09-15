"""A separately charged, post-freeze same-family magnetic continuation.

This diagnostic changes no operators or exact constraints. It resumes a saved
optimization in a new case, preserving the original frozen attempt and clock.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import time
from research.interacting_scaling_20260915.budget import ROOT, OUT, STD, NUM, dump
from research.interacting_scaling_20260915.pipeline import execute, PREFIX


def initialize(source, name):
    frame = json.loads((source/'prepared/frame.json').read_text())
    design = json.loads((source/'design.json').read_text())
    if frame.get('magnetization') != 1 or design.get('magnetization') != 1:
        raise ValueError('Only continuation within the same M_S=1 family is supported')
    if frame['fixture_sha256'] != hashlib.sha256((source/'fixture.json').read_bytes()).hexdigest():
        raise ValueError('Prepared maps belong to a different Hamiltonian')
    checkpoint = source/'solve/checkpoint.npz'
    if not checkpoint.exists() or not (source/'exact.json').exists():
        raise ValueError('A completed saved optimization with an exact screen check is required')
    target = OUT/'cases'/name
    target.mkdir(parents=True, exist_ok=False)
    for item in ('fixture.json', 'upper.json', 'design.json'):
        shutil.copyfile(source/item, target/item)
    shutil.copytree(source/'prepared', target/'prepared')
    shutil.copyfile(checkpoint, target/'resume.npz')
    files = [source/'fixture.json', source/'design.json', checkpoint]+sorted((source/'prepared').iterdir())
    dump(target/'continuation_dependency.json', {'source_case': str(source),
        'all_source_discovery_and_failed_attempt_costs_additional': True,
        'post_freeze_optimization_diagnostic_not_a_frozen_heldout_result': True,
        'same_Hamiltonian_sector_operators_constraints_and_maps': True,
        'saved_primal_and_dual_iterates_preserved': True,
        'source_sha256': {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in files if p.is_file()}})
    return target


def run(source, name, seconds=300, deadline=None):
    if type(seconds) is not int or not 1 <= seconds <= 840:
        raise ValueError('A bounded continuation of at most 840 solver seconds is required')
    last = json.loads((source/'solve/history.json').read_text())[-1]
    mu = float(last['mu'])
    if mu not in (2., .03): raise ValueError('Unsupported saved refinement phase')
    start = time.monotonic()
    target = OUT/'cases'/name
    result = execute(name, [('initialize', 30, [STD, '-B', '-S', '-m', PREFIX+'screen_continue', 'initialize', str(source), name]),
        ('solve', seconds+45, [NUM, '-B', '-m', PREFIX+'solve', str(target), 'solve', '--seconds', str(seconds), '--mu', str(mu), '--restart', str(target/'resume.npz')]),
        ('check', 180, [STD, '-B', '-S', '-m', PREFIX+'nonsinglet', 'check', str(target), 'solve', str(target/'exact.json')])], deadline)
    record = {'source_case': str(source), 'execution': result,
        'additional_continuation_seconds': time.monotonic()-start,
        'fresh_complete_cold_calculation': False, 'earlier_costs_additional': True}
    if result['completed']:
        record['exact_screen'] = json.loads((target/'exact.json').read_text())
    dump(target/'continuation_result.json', record)
    return record


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('action', choices=('initialize', 'run'))
    p.add_argument('source', type=Path); p.add_argument('name'); p.add_argument('--seconds', type=int, default=300)
    a = p.parse_args()
    initialize(a.source.resolve(), a.name) if a.action == 'initialize' else run(a.source.resolve(), a.name, a.seconds)
