"""Complete comparison export repairs and the predeclared initialization control."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import time
from research.transfer_solver_20260915.budget import ROOT, OUT, dump
from research.transfer_solver_20260915.campaign import STD, NUM, CHEM, PREFIX


def bounded(name, seconds, command):
    return subprocess.run([STD, '-B', '-S', '-m', PREFIX+'budget', '--name', name,
                            '--seconds', str(seconds), '--']+command, cwd=ROOT).returncode


def checkpoint_diagnostic(case):
    from experiments.marginal_symbolic import decode
    from research.reconstruction_compression_20260914.moments import Proposal
    start = time.monotonic()
    data = json.loads((case/'fixture.json').read_text())
    path = case/'mps/checkpoint.json'
    state = json.loads(path.read_text())
    h = decode(data['hamiltonian'], data['modes'], 4)
    proposal = Proposal(data, state)
    proposal.fill(h)
    e = sum(float(c)*proposal.cache[w] for w, c in h.items())
    reference = json.loads((case/'fci_reference.json').read_text())
    dump(case/'checkpoint_energy_proposal.json', {'numerical_projected_checkpoint_energy_Ha': e,
        'numerical_reference_difference_mHa': 1000*(e-reference['numerical_electronic_energy_Ha']),
        'numerical_projected_norm': proposal.norm, 'checkpoint_sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
        'seconds': time.monotonic()-start, 'certified_upper': False, 'used_for_discovery': False,
        'purpose': 'Diagnose the actual charge-projected checkpoint, separately from penalized DMRG energy'})


def fresh_benchmark(case):
    from research.transfer_solver_20260915 import benchmark
    from research.transfer_followup_20260915.repair_comparisons import typed_dump
    benchmark.dump = typed_dump
    benchmark.run(case)


def complete_checkpoint(case):
    execution = json.loads((case/'execution.json').read_text())
    if execution['status'] == 'completed' or execution['last_step'] != 'state':
        return
    checkpoint = case/'mps/checkpoint.json'
    if not checkpoint.exists():
        return
    protocol = json.loads((case/'protocol.json').read_text())
    destination = case/'mps/state.json'
    if destination.exists():
        raise ValueError('Refusing to replace an exported state')
    shutil.copyfile(checkpoint, destination)
    steps = [row for row in protocol['steps'] if row['name'] not in ('generate', 'state')]
    dump(case/'checkpoint_completion_protocol.json', {'checkpoint_sha256': hashlib.sha256(checkpoint.read_bytes()).hexdigest(),
        'source_is_this_fresh_runs_completed_sweep': True, 'DMRG_convergence_claimed': False,
        'exact_upper_replay_required': True, 'original_failed_state_job_remains_charged': True, 'steps': steps,
        'reason': 'Any structurally valid nonzero fixed-N state can supply an upper; convergence is not required for Rayleigh certification.'})
    started = time.monotonic()
    status = 'completed'
    for row in steps:
        code = bounded('h10_hf64_checkpoint_'+row['name'], row['seconds'], row['command'])
        if code:
            status = 'stopped_after_failed_step'
            break
    dump(case/'checkpoint_completion.json', {'status': status, 'last_step': row['name'],
        'seconds_after_saved_sweep': time.monotonic()-started, 'failed_discovery_cost_is_additional': True})


def run():
    if not (OUT/'followup_execution.json').exists():
        raise ValueError('Run these refinements after the declared comparison set')
    os.environ['NUMBA_CACHE_DIR'] = str(OUT/'refinement_numba_cache')
    results = []
    repair = subprocess.run([STD, '-B', '-S', '-m', 'research.transfer_followup_20260915.repair_comparisons', 'run'], cwd=ROOT)
    results.append({'name': 'FCI_scalar_export_repair', 'exit_code': repair.returncode})
    frozen = json.loads((OUT/'hf_size_predeclared.json').read_text())
    if hashlib.sha256((ROOT/frozen['source']).read_bytes()).hexdigest() != frozen['source_sha256']:
        raise ValueError('HF initialization rule changed after preregistration')
    source = OUT/'adaptive/h10_bond64'
    if (source/'execution.json').exists() and json.loads((source/'execution.json').read_text())['status'] != 'completed':
        result = subprocess.run([STD, '-B', '-S', '-m', 'research.transfer_followup_20260915.hf_size', 'run'], cwd=ROOT)
        results.append({'name': 'HF_product_initialization_control', 'exit_code': result.returncode})
        case = OUT/'adaptive/h10_hf64'
        complete_checkpoint(case)
        if (case/'fixture.json').exists():
            results.append({'name': 'h10_hf64_fci', 'exit_code': bounded('h10_hf64_fci', 180,
                [CHEM, '-B', '-m', 'research.transfer_followup_20260915.refinements', 'benchmark', str(case)])})
    for case in (OUT/'cases/h10_size', OUT/'adaptive/h10_bond64', OUT/'adaptive/h10_hf64'):
        if (case/'mps/checkpoint.json').exists() and (case/'fci_reference.json').exists():
            results.append({'name': case.name+'_checkpoint_diagnostic', 'exit_code': bounded(case.name+'_checkpoint_diagnostic', 90,
                [NUM, '-B', '-m', 'research.transfer_followup_20260915.refinements', 'diagnose', str(case)])})
    dump(OUT/'refinement_execution.json', results)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('action', choices=('run', 'diagnose', 'benchmark'))
    p.add_argument('case', type=Path, nargs='?')
    a = p.parse_args()
    if a.action == 'run':
        run()
    elif a.action == 'diagnose':
        checkpoint_diagnostic(a.case.resolve())
    else:
        fresh_benchmark(a.case.resolve())
