"""Transport low-order MPS moments to the original orbital basis for discovery.

Only reduced operator arrays through degree six are transformed. The numerical
moments guide optimization; the unchanged exact lower checker accepts the proof.
"""
import argparse
from fractions import Fraction as F
from itertools import combinations
import hashlib
import json
from math import comb
from pathlib import Path
import shutil
import subprocess
import time
from research.transfer_solver_20260915.budget import ROOT, OUT, dump
from research.transfer_solver_20260915.campaign import STD, NUM, PREFIX

SOURCE = OUT/'adaptive/h10_rotated_upper'
LOCAL = OUT/'rotated_h10'
CASE = OUT/'adaptive/h10_correlated_guide'


def compound(O, sets):
    import numpy as np
    indices = np.asarray(sets, dtype=int)
    return np.linalg.det(O[indices[:, None, :, None], indices[None, :, None, :]])


def canonical_moment(word, matrices, indices):
    if not word:
        return 1.
    left = tuple(i for c, i in word if c)
    right = tuple(i for c, i in word if not c)
    if len(left) != len(right) or sum(i % 2 == 0 for i in left) != sum(i % 2 == 0 for i in right):
        return 0.
    degree = len(left)
    key = degree, sum(i % 2 == 0 for i in left)
    sign = (-1)**(degree*(degree-1)//2)
    return sign*float(matrices[key][indices[key][left], indices[key][right]])


def prepare():
    import numpy as np
    from scipy import sparse
    from scipy.sparse.linalg import lsqr
    from experiments.marginal_symbolic import decode
    from research.reconstruction_compression_20260914.moments import Proposal
    from research.molecular_collective_20260913.core import digest
    start = time.monotonic()
    read = lambda p: json.loads(p.read_text())
    data, localized = read(SOURCE/'fixture.json'), read(LOCAL/'fixture.json')
    rotation = read(LOCAL/'rotation.json')
    if rotation['original_fixture_sha256'] != digest(data):
        raise ValueError('Moment transformation changed the original model')
    state = read(LOCAL/'mps/state.json')
    oracle = Proposal(localized, state)
    den = int(rotation['denominator'])
    V = np.array([[float(F(v, den)) for v in row] for row in rotation['integer_matrix']])
    O = np.kron(V, np.eye(2))
    matrices, indices, checks = {}, {}, []
    for degree in range(1, 4):
        for alpha in range(degree+1):
            sets = [x for x in combinations(range(data['modes']), degree) if sum(i % 2 == 0 for i in x) == alpha]
            key = degree, alpha
            words = {(i, j): tuple((1, p) for p in a)+tuple((0, q) for q in reversed(b))
                     for i, a in enumerate(sets) for j, b in enumerate(sets[:i+1])}
            oracle.fill(words.values())
            local = np.zeros((len(sets), len(sets)))
            for (i, j), w in words.items():
                local[i, j] = local[j, i] = oracle.cache[w]
            W = compound(O, sets)
            canonical = W@local@W.T
            expected = comb(state['spin_counts'][0], alpha)*comb(state['spin_counts'][1], degree-alpha)
            trace_error = abs(float(np.trace(canonical))-expected)
            if trace_error > 1e-7*max(1, expected):
                raise ValueError('Transformed RDM trace violates fixed charge')
            matrices[key], indices[key] = canonical, {v: i for i, v in enumerate(sets)}
            checks.append({'degree': degree, 'alpha_count': alpha, 'dimension': len(sets),
                'trace_error': trace_error, 'seconds': time.monotonic()-start})
            print(json.dumps(checks[-1]), flush=True)
    h = decode(data['hamiltonian'], data['modes'], 4)
    energy = sum(float(v)*canonical_moment(w, matrices, indices) for w, v in h.items())
    upper = read(LOCAL/'original_upper.json')
    gap = abs(energy-float(F(upper['upper_Ha'])))
    if gap > 1e-6:
        raise ValueError(('Transported physical energy disagrees with the exact upper', gap))
    prepared = CASE/'prepared'
    shutil.copytree(SOURCE/'prepared', prepared)
    meta = read(prepared/'frame.json')
    rows = [tuple(tuple(v) for v in w) for w in meta['rows']]
    weights = np.load(prepared/'weights.npy')
    physical = weights*np.array([canonical_moment(w, matrices, indices) for w in rows])
    T = sparse.load_npz(prepared/'twirl.npz')
    selected = np.load(prepared/'selected.npy')
    pullback = lsqr(T[selected].T, T.T@physical, atol=1e-13, btol=1e-13, iter_lim=1000)
    np.save(prepared/'physical_dual.npy', pullback[0]/np.load(prepared/'scale.npy'))
    np.save(prepared/'moments.npy', physical)
    meta['physical_dual_override'] = 'Degree-1/2/3 RDMs of this pass local MPS transported to the original basis'
    meta['physical_dual_source_MPS_sha256'] = hashlib.sha256((LOCAL/'mps/state.json').read_bytes()).hexdigest()
    (prepared/'frame.json').write_text(json.dumps(meta, indent=2)+'\n')
    dump(CASE/'moment_transport.json', {'source_MPS': str(LOCAL/'mps/state.json'),
        'source_prepared': str(SOURCE/'prepared'), 'numerical_proposal_only': True,
        'full_N_determinants_enumerated': 0, 'highest_reduced_density_order': 3,
        'checks': checks, 'transported_energy_Ha': energy, 'absolute_difference_from_exact_upper_Ha': gap,
        'moment_pullback_residual': float(pullback[3]), 'MPS_transfer_steps': oracle.steps,
        'local_words_evaluated': len(oracle.cache), 'seconds': time.monotonic()-start,
        'old_maps_shared_cost_additional': True, 'initial_Gram_matrices': 'zero'})


def run():
    CASE.mkdir(exist_ok=False)
    for name in ('fixture.json', 'upper.json', 'nonsinglet.json'):
        shutil.copyfile(SOURCE/name, CASE/name)
    steps = [('prepare', 300, [NUM, '-B', '-m', 'research.transfer_followup_20260915.correlated_guide', 'prepare']),
        ('stage1', 900, [NUM, '-B', '-m', PREFIX+'solve', str(CASE), 'stage1', '--seconds', '780', '--mu', '2']),
        ('stage2', 450, [NUM, '-B', '-m', PREFIX+'solve', str(CASE), 'stage2', '--seconds', '350', '--mu', '.03',
                       '--restart', str(CASE/'stage1/checkpoint.npz')]),
        ('replay', 600, [STD, '-B', '-S', '-m', PREFIX+'actions', str(CASE), 'replay'])]
    dump(CASE/'protocol.json', {'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'reason': 'Large HF-guided residuals and measured 161-second setup cost in preceding H10 solve',
        'new_physical_dual_from_this_pass_correlated_state': True, 'larger_recorded_solve_budget': True,
        'same_family_and_exact_checker': True, 'initial_Gram_and_ideal_coordinates': 'zero',
        'previous_maps_upper_and_nonsinglet_costs_additional': True, 'FCI_input': False, 'steps': steps})
    started = time.monotonic()
    outcomes = []
    for name, seconds, command in steps:
        result = subprocess.run([STD, '-B', '-S', '-m', PREFIX+'budget', '--name', 'h10_correlated_'+name,
            '--seconds', str(seconds), '--']+command, cwd=ROOT)
        outcomes.append({'stage': name, 'exit_code': result.returncode})
        if result.returncode:
            break
    dump(CASE/'execution.json', {'steps': outcomes, 'elapsed_seconds': time.monotonic()-started,
        'last_step': name, 'status': 'completed' if len(outcomes) == len(steps) and not outcomes[-1]['exit_code'] else 'stopped_after_failed_step'})


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('action', choices=('run', 'prepare'))
    a = p.parse_args()
    run() if a.action == 'run' else prepare()
