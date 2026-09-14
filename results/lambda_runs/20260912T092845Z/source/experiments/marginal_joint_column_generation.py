"""Restricted-master discovery in a generated coefficient positivity dictionary.

The full cached dictionary is currently retained for numerical pricing. This
reduces LP working columns, not dictionary generation or quotient dimension.
Only the existing full rational export establishes certificate acceptance.
"""
from fractions import Fraction as F
from pathlib import Path
import json
import time

from experiments.marginal_joint_coefficient_constructor import input_digest, export


def construct(data, prepared, out, batch=32, max_rounds=80, time_limit=120,
              metric_floor=F(9, 10000), numerator_floor=F(1, 10000)):
    import numpy as np
    from scipy.sparse import load_npz, hstack, eye
    from scipy.optimize import linprog
    if type(batch) is not int or batch < 1 or type(max_rounds) is not int or max_rounds < 1 or time_limit <= 0:
        raise ValueError('Positive batch, rounds and time budget required')
    prepared, out = Path(prepared), Path(out)
    meta = json.loads((prepared/'prepared.json').read_text())
    if set(data) != {'modes', 'particles', 'hamiltonian'} or meta['source_sha256'] != input_digest(data):
        raise ValueError('Prepared dictionary is not bound to bare Hamiltonian')
    matrix = load_npz(prepared/'model.npz').tocsc()
    rhs = np.load(prepared/'rhs.npy', allow_pickle=False)
    n, a = len(meta['orbits']), len(meta['atom_orbits'])
    rows = len(rhs)
    if matrix.shape != (rows, n+2*a+1) or rows % 2 != 1:
        raise ValueError('Malformed joint dictionary')
    rhs[:(rows-1)//2] *= float(metric_floor)/.001
    rhs -= matrix[:, -1].toarray()[:, 0]*float(numerator_floor)
    metric = matrix[:, :n]
    atoms = matrix[:, n:-1]
    # Geometry-only seeds: no solution, metric coefficients or active atom list.
    seed = [j for j, orbit in enumerate(meta['atom_orbits'])
            if min(meta['labels'][member][1].bit_count() for member in orbit) <= 2]
    active = sorted(seed + [a+j for j in seed])
    inactive = np.ones(2*a, dtype=bool); inactive[active] = False
    identity = eye(rows, format='csc')
    history = []; started = time.monotonic(); accepted = False
    out.mkdir(parents=True, exist_ok=True)
    reason = 'Round budget reached'
    for iteration in range(max_rounds):
        remaining = time_limit-(time.monotonic()-started)
        if remaining <= 0:
            reason = 'Total time budget reached'; break
        master = hstack((metric, atoms[:, active], identity, -identity), format='csc')
        objective = np.r_[np.zeros(n+len(active)), np.ones(2*rows)]
        result = linprog(objective, A_eq=master, b_eq=rhs,
                         bounds=[(None, None)]*n+[(0, None)]*(len(active)+2*rows),
                         method='highs-ds', options={'time_limit': remaining,
                         'dual_feasibility_tolerance': 1e-8, 'primal_feasibility_tolerance': 1e-8})
        record = {'round': iteration, 'active_atoms': len(active), 'master_columns': master.shape[1],
                  'master_nonzeros': master.nnz, 'solver_status': result.message,
                  'elapsed_seconds': time.monotonic()-started}
        if not result.success:
            reason = 'Restricted solve did not finish'; history.append(record); break
        x = result.x
        record['phase_one_l1'] = float(result.fun)
        # For zero-cost nonnegative columns, reduced cost = -A_j^T y.
        scores = np.asarray(atoms.T @ result.eqlin.marginals)
        record['maximum_inactive_dual_violation'] = float(max(scores[inactive], default=0.))
        if result.fun < 1e-7:
            proposal = {'source_sha256': input_digest(data), 'orbits': meta['orbits'],
                        'metric_coefficients': x[:n].tolist(), 'metric_bound': float(metric_floor),
                        'numerator_bound': float(numerator_floor), 'gamma': meta['gamma']}
            weights = np.zeros(2*a); weights[active] = x[n:n+len(active)]
            for name, values in [('weight', weights[:a]), ('numerator', weights[a:])]:
                pairs = [(meta['labels'][member], float(value)/len(orbit))
                         for value, orbit in zip(values, meta['atom_orbits']) if abs(value) > 1e-13
                         for member in orbit]
                proposal[name+'_labels'] = [p[0] for p in pairs]
                proposal[name+'_values'] = [p[1] for p in pairs]
            (out/'proposal.json').write_text(json.dumps(proposal, indent=2)+'\n')
            try:
                proof = export(data, proposal, out/'proof')
                accepted = True; reason = 'Exact certificate accepted'
                record['exact_receipt'] = proof
            except ValueError as error:
                record['exact_export_rejection'] = str(error)
                reason = 'Numerically feasible proposal failed exact export'
            history.append(record)
            print(json.dumps(record), flush=True)
            break
        selected = []
        for block in (0, a):
            candidates = np.flatnonzero(inactive[block:block+a] & (scores[block:block+a] > 1e-8))+block
            order = sorted(candidates, key=lambda j: (-scores[j], j))[:batch]
            selected.extend(map(int, order))
        record['added_atoms'] = len(selected)
        history.append(record)
        print(json.dumps(record), flush=True)
        (out/'history.json').write_text(json.dumps(history, indent=2)+'\n')
        if not selected:
            reason = 'No numerically improving dictionary column; no exact infeasibility claim'; break
        active.extend(selected); inactive[selected] = False
    receipt = {'source_sha256': input_digest(data), 'exact_accepted': accepted, 'reason': reason,
               'elapsed_seconds_including_export': time.monotonic()-started,
               'dictionary_columns': 2*a, 'metric_columns': n, 'coefficient_rows': rows,
               'initial_atom_columns': 2*len(seed), 'final_active_atom_columns': len(active),
               'rounds': len(history), 'batch_per_positivity_block': batch,
               'scope': 'Generated geometry-only seed and dual pricing over complete cached dictionary. No inherited solution or selected proof atoms. Artificial coefficient residuals are proposal-only; full exact export required. Complete dictionary storage and quotient dimension remain.'}
    (out/'receipt.json').write_text(json.dumps(receipt, indent=2)+'\n')
    (out/'history.json').write_text(json.dumps(history, indent=2)+'\n')
    (out/'active_columns.json').write_text(json.dumps(active)+'\n')
    print(json.dumps(receipt), flush=True)
    return receipt


if __name__ == '__main__':
    import argparse
    from unittest.mock import patch
    parser = argparse.ArgumentParser()
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--prepared', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--batch', type=int, default=32)
    parser.add_argument('--time-limit', type=float, default=120)
    args = parser.parse_args()
    with patch('experiments.marginal_determinant_tree.DeterminantOracle.action', side_effect=AssertionError('No states')), \
         patch('experiments.marginal_spin_constructor.spin_states', side_effect=AssertionError('No state generation')), \
         patch('experiments.marginal_polynomial_metric.complete_number_ideals', side_effect=AssertionError('No full-population lift')):
        construct(json.loads(args.source.read_text()), args.prepared, args.out, batch=args.batch, time_limit=args.time_limit)
