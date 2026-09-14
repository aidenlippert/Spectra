"""Bounded dual-mixture discovery for two added coefficients with old recipe frozen."""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import importlib
import json
import math
import sys
import time
import numpy as np
from scipy.linalg import eigh
from scipy.optimize import minimize, linprog
from spin_word_fraction_free import solve

from residual_coherence_obstruction import ROOT, BASE, OLD_MODULES, pure_coherence
from spin_word_numeric import prepare, legacy_terms
from experiments.marginal_spin_word_telescope import local_value
from experiments.marginal_hopping_telescope import projected_matrix


def action_from_matrix(matrix, denominator):
    result = [{} for _ in range(4096)]
    for (r, c), value in matrix.items():
        result[c][r] = F(value, denominator)
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('seed', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    certificate = json.loads(args.seed.read_text())
    candidates_path = BASE.parent / 'joint_projector/signed_density/symmetry_diagonal_candidates.json'
    candidates = json.loads(candidates_path.read_text())['candidates']
    shapes = [{int(s): v for s, v in candidates[i]['diagonal'].items()} for i in (1, 2, 4, 5, 6, 7, 31, 11, 20)]
    x, physical, blocks, *_ = prepare(certificate, shapes)
    old_actions, coefficients = [], []
    for name in OLD_MODULES:
        module = importlib.import_module('experiments.marginal_' + name)
        field = 'coherent_projector' if name == 'coherent_projector_telescope' else name
        for label in getattr(module, 'LABELS', ('single',)):
            old_actions.append(module.actions({label: 1}) if hasattr(module, 'LABELS') else module.actions())
            coefficients.append(F(certificate[field].get(label, 0)) if hasattr(module, 'LABELS') else F(certificate[field]))
    terms = legacy_terms(certificate)
    for label, value in certificate['spin_word_telescope'].items():
        terms[label] += F(value)
    alpha, beta = F(certificate['penalty']), F(certificate['joint']['penalty'])
    theta_h = F(certificate['projector_sum_ceiling']) / certificate['windows']
    theta_j = F(certificate['joint']['projector_sum_ceiling']) / certificate['joint']['windows']
    penalty_cost = alpha * theta_h + beta * theta_j
    new_actions, receipts = [], []
    for case in ('W_zero', 'W_plus_1'):
        path = BASE / case / 'telescope_replay.json'
        receipt = json.loads(path.read_text())
        if receipt.get('accepted') is not True:
            raise ValueError('Accepted new telescope required')
        for name, digest in receipt['source_sha256'].items():
            if hashlib.sha256((ROOT / name).read_bytes()).hexdigest() != digest:
                raise ValueError('Stale telescope source')
        vector = {int(s): a for s, a in receipt['five_site_vector'].items()}
        denominator, _, matrix, _ = pure_coherence(vector)
        new_actions.append(action_from_matrix(matrix, denominator))
        receipts.append(path)
    fresh = physical(x)
    matrices = []
    discrepancy = 0.0
    for key, a, ds, ph, q, ts, cols, *_ in blocks:
        scale = np.sqrt([sum(v*v for v in column.values()) for column in cols])
        old = np.array([projected_matrix(cols, action) for action in old_actions], float) / scale[:, None] / scale[None, :]
        diagonal = np.array([float(local_value(next(iter(column)), terms)) for column in cols])
        base = a + np.einsum('i,ijk->jk', np.array(coefficients, float), old) + np.diag(diagonal) + float(alpha)*ph + float(beta)*q
        direct = np.array(fresh[key][1], float) / scale[:, None] / scale[None, :]
        direct += np.einsum('i,ijk->jk', np.array(coefficients, float), old) + np.diag(diagonal) + float(alpha)*ph + float(beta)*q
        discrepancy = max(discrepancy, float(np.max(abs(direct - base))))
        new = np.array([projected_matrix(cols, action) for action in new_actions], float) / scale[:, None] / scale[None, :]
        matrices.append((key, base, new))
    if discrepancy > 1e-9:
        raise ValueError('Fresh physical matrix reconstruction disagrees')
    columns = {key: cols for key, a, ds, ph, q, ts, cols, *_ in blocks}
    rows, energies, states = [], [], []
    seen = set()

    def add_atom(key, vector, base):
        cols = columns[key]
        scale = np.sqrt([sum(v*v for v in col.values()) for col in cols])
        unnormalized = vector / scale
        coeff = [round(float(v/max(abs(unnormalized)))*10**7) for v in unnormalized]
        state = {s: z*a for z, col in zip(coeff, cols) if z for s, a in col.items()}
        canonical = tuple(sorted(state.items()))
        if canonical in seen:
            return
        seen.add(canonical)
        norm = sum(a*a for a in state.values())
        exact = [F(1)]
        for action in new_actions:
            exact.append(sum((F(a)*value*state.get(target, 0) for s, a in state.items()
                              for target, value in action[s].items()), F(0))/norm)
        w = np.array(coeff)*scale/np.sqrt(float(norm))
        rows.append(exact)
        energies.append(float(w @ base @ w))
        states.append({str(s): a for s, a in state.items()})

    for key, base, new in matrices:
        ev, vectors = eigh(base, subset_by_index=[0, 0])
        add_atom(key, vectors[:, 0], base)
    pricing = []
    for iteration in range(21):
        scaling = np.array([1., 1000., 1000.])
        array = np.array(rows, float).T
        result = linprog(energies, A_eq=scaling[:, None]*array, b_eq=[1., 0., 0.], bounds=(0, None),
                         method='highs-ds', options={'presolve': False,
                         'dual_feasibility_tolerance': 1e-10, 'primal_feasibility_tolerance': 1e-10})
        if not result.success:
            raise ValueError('Numerical three-row family LP failed: ' + result.message)
        if iteration == 20:
            break
        dual = result.eqlin.marginals*scaling
        lowest = float('inf')
        count = len(states)
        for key, base, new in matrices:
            ev, vectors = eigh(base - np.einsum('i,ijk->jk', dual[1:], new), subset_by_index=[0, 0])
            reduced = float(ev[0]) - dual[0]
            lowest = min(lowest, reduced)
            if reduced < -1e-12:
                add_atom(key, vectors[:, 0], base)
        pricing.append({'round': iteration, 'candidate_upper': (float(result.fun)-float(penalty_cost))/5,
                        'minimum_reduced_eigenvalue': lowest, 'added': len(states)-count,
                        'dual_new_coefficients': list(map(float, -dual[1:]))})
        if len(states) > 2500:
            raise ValueError('Candidate bound exceeded')
        if count == len(states):
            break
    support = [i for i, weight in enumerate(result.x) if weight > 1e-14]
    output = {'accepted': False, 'proposal_written': False, 'candidate_count': len(states),
              'basis_size': len(support), 'pricing': pricing,
              'numerical_upper': (float(result.fun)-float(penalty_cost))/5,
              'seconds': time.monotonic()-started,
              'scope': 'Numerical selection only. At most3 exact positive sources cancelling both new moments can cap arbitrary new coefficients with all preceding fields frozen. Independent physical replay is mandatory.'}
    try:
        weights = solve([[rows[i][j] for i in support] for j in range(3)], [1, 0, 0])
        if min(weights) < 0:
            raise ValueError('Negative exact proposed weight')
        proposal = {'kind': 'residual_coherence_fixed_family_proposal_v1',
                    'seed_certificate': str(args.seed.resolve().relative_to(ROOT)),
                    'mixture': [{'weight': str(weight), 'vector': states[i]} for i, weight in zip(support, weights) if weight],
                    'scope': 'Untrusted mixture proposed to cap only the two new coefficients; all old fields fixed.'}
        (args.output / 'family_proposal.json').write_text(json.dumps(proposal, indent=2) + '\n')
        output['proposal_written'] = True
    except ValueError as exc:
        output['exact_failure'] = str(exc)
    files = {Path(__file__).resolve(), args.seed.resolve(), candidates_path, *receipts}
    for module in tuple(sys.modules.values()):
        name = getattr(module, '__file__', None)
        if name:
            path = Path(name).resolve()
            if path.is_relative_to(ROOT / 'experiments') or path.is_relative_to(BASE.parent / 'discovery'):
                files.add(path)
    output['source_sha256'] = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(files)}
    (args.output / 'family_diagnostic.json').write_text(json.dumps(output, indent=2) + '\n')
    print(json.dumps({k: v for k, v in output.items() if k not in ('source_sha256', 'pricing', 'scope')}))


if __name__ == '__main__':
    main()
