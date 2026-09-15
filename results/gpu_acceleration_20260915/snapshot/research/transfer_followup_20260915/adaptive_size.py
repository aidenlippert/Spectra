"""Declared smaller-state retry only after the frozen size run hits its state budget."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import time
from research.transfer_solver_20260915.budget import ROOT, OUT, dump
from research.transfer_solver_20260915.campaign import STD, NUM, CHEM, PREFIX, CASES


def state(case):
    from research.correlated_pair_20260913.mps_spatial import run
    data = json.loads((case/'fixture.json').read_text())
    output = case/'mps'
    output.mkdir(exist_ok=False)
    run(data, [data['particles']//2]*2, output, bond=64, sweeps=6, initial=None)


def prepare(case):
    from research.collective_completion_20260914 import spin_rows
    from research.sector_quotient_20260914.fast_twirl import twirl
    from research.transfer_solver_20260915.prepare import build
    # Proposal-side replacement; the original independent exact accepting path remains unchanged.
    spin_rows.twirl = twirl
    build(case)


def run():
    execution = json.loads((OUT/'cases/h10_size/execution.json').read_text())
    if execution['status'] == 'completed' or execution['last_step'] not in ('state', 'upper'):
        raise ValueError('This retry specifically addresses a measured state/upper construction failure')
    case = OUT/'adaptive/h10_bond64'
    case.mkdir(parents=True, exist_ok=False)
    dump(case/'specification.json', {**CASES['h10_size'], 'role': 'exploratory retry after frozen state-budget failure'})
    steps = [('generate', 90, [CHEM, '-B', '-m', PREFIX+'generate', str(case)]),
             ('state', 600, [NUM, '-B', '-m', 'research.transfer_followup_20260915.adaptive_size', 'state', str(case)]),
             ('upper', 300, [STD, '-B', '-S', '-m', PREFIX+'actions', str(case), 'upper']),
             ('nonsinglet', 240, [NUM, '-B', '-m', PREFIX+'actions', str(case), 'nonsinglet']),
             ('prepare', 900, [NUM, '-B', '-m', 'research.transfer_followup_20260915.adaptive_size', 'prepare', str(case)]),
             ('stage1', 400, [NUM, '-B', '-m', PREFIX+'solve', str(case), 'stage1', '--seconds', '300', '--mu', '2']),
             ('stage2', 300, [NUM, '-B', '-m', PREFIX+'solve', str(case), 'stage2', '--seconds', '180', '--mu', '.03',
                              '--restart', str(case/'stage1/checkpoint.npz')]),
             ('replay', 600, [STD, '-B', '-S', '-m', PREFIX+'actions', str(case), 'replay'])]
    sources = json.loads((OUT/'validation_set_protocol.json').read_text())['source_hashes']
    sources[str(Path(__file__).relative_to(ROOT))] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    dump(case/'protocol.json', {'original_frozen_failure': execution, 'fresh_integrals_and_random_state': True,
        'bond_cap': 64, 'sweeps_cap': 6, 'previous_checkpoint_used': False,
        'proposal_projector': 'Existing exact combinatorial spin twirl', 'accepting_checker_unchanged': True,
        'same_lower_dictionary_and_optimizer_as_frozen_campaign': True, 'frozen_source_hashes': sources,
        'two_declared_changes': ['smaller MPS bond/sweep budget', 'faster proposal-side exact twirl'],
        'not_an_untouched_frozen_validation_case': True,
        'steps': [{'name': n, 'seconds': t, 'command': c} for n, t, c in steps]})
    start = time.monotonic()
    status = 'completed'
    for name, seconds, command in steps:
        if any(hashlib.sha256((ROOT/p).read_bytes()).hexdigest() != digest for p, digest in sources.items()):
            raise ValueError('Adaptive-case source changed while running')
        result = subprocess.run([STD, '-B', '-S', '-m', PREFIX+'budget', '--name', 'h10_bond64_'+name,
                                 '--seconds', str(seconds), '--']+command, cwd=ROOT)
        if result.returncode:
            status = 'stopped_after_failed_step'
            break
    dump(case/'execution.json', {'status': status, 'last_step': name, 'elapsed_seconds': time.monotonic()-start,
        'original_failed_run_cost_is_additional': True, 'new_from_integrals_attempt': True})


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('action', choices=('run', 'state', 'prepare'))
    p.add_argument('case', type=Path, nargs='?')
    a = p.parse_args()
    if a.action == 'run':
        run()
    elif a.action == 'state':
        state(a.case.resolve())
    else:
        prepare(a.case.resolve())
