"""Spin-exchange reduction of an untrusted bare-H coefficient LP proposal.

Exact compiled-polynomial symmetry is checked before numerical projection.
The original unsymmetrized rational exporter remains the acceptance gate.
"""
from fractions import Fraction as F
from pathlib import Path
import json
import time

from experiments.marginal_joint_coefficient_constructor import (
    CoefficientQuotient, input_digest, vanishing_feature, export,
)
from experiments.marginal_polynomial_metric import JointPolynomial


def flip_mask(mask, sites):
    if type(sites) is not int or sites < 1 or type(mask) is not int or not 0 <= mask < 1 << (2*sites):
        raise ValueError('Valid occupation mask and positive site count required')
    even = sum(1 << (2*i) for i in range(sites))
    return ((mask & even) << 1) | ((mask & (even << 1)) >> 1)


def flip_label(label, sites):
    family, required, occupied = label
    return [family, flip_mask(required, sites), flip_mask(occupied, sites)]


def _orbits(permutation):
    size = len(permutation)
    if sorted(permutation) != list(range(size)) or any(permutation[permutation[i]] != i for i in range(size)):
        raise ValueError('An involutive permutation is required')
    return [[i] if j == i else [i, j] for i, j in enumerate(permutation) if i <= j]


def _average_projection(orbits, size):
    from scipy.sparse import csc_matrix
    ri, ci, values = [], [], []
    for i, orbit in enumerate(orbits):
        for member in orbit:
            ri.append(i); ci.append(member); values.append(1/len(orbit))
    return csc_matrix((values, (ri, ci)), shape=(len(orbits), size))


def reduce_cached(model_dir, source, out):
    import numpy as np
    from scipy.sparse import load_npz, save_npz, block_diag
    model_dir, out = Path(model_dir), Path(out)
    data = json.loads(Path(source).read_text())
    meta = json.loads((model_dir/'prepared.json').read_text())
    if set(data) != {'modes', 'particles', 'hamiltonian'} or meta.get('source_sha256') != input_digest(data):
        raise ValueError('Cached model is not bound to bare source')
    matrix = load_npz(model_dir/'model.npz').tocsc()
    rhs = np.load(model_dir/'rhs.npy', allow_pickle=False)
    orbits, labels = meta['orbits'], meta['labels']
    n, a, sites = len(orbits), len(labels), data['modes']//2
    ring = JointPolynomial(dict(data, polynomial_metric=vanishing_feature([[0]*sites], sites)))
    quotient = CoefficientQuotient(ring, 6)
    qsize = len(quotient.basis)
    if matrix.shape != (2*qsize+1, n+2*a+1) or rhs.shape != (2*qsize+1,):
        raise ValueError('Malformed cached joint model')
    qp = [quotient.index[flip_mask(m, sites)] for m in quotient.basis]
    row_orbits = _orbits(qp)
    lookup = {tuple(label): i for i, label in enumerate(labels)}
    if len(lookup) != a:
        raise ValueError('Duplicate atom labels')
    try:
        ap = [lookup[tuple(flip_label(label, sites))] for label in labels]
    except KeyError as error:
        raise ValueError('Atom dictionary is not spin-flip closed') from error
    atom_orbits = _orbits(ap)
    atoms = matrix[:qsize, n:n+a].tocsc()
    if (atoms[qp, :][:, ap] != atoms).nnz:
        raise ValueError('Cached atoms violate coefficient spin-flip equivariance')
    if (matrix[qsize:2*qsize, n+a:n+2*a] != atoms).nnz:
        raise ValueError('Cached atom blocks disagree')
    # Compile all generated metric directions exactly; no physical assignments.
    started = time.monotonic()
    for j, orbit in enumerate(orbits):
        trial = JointPolynomial(dict(data, polynomial_metric=vanishing_feature(orbit, sites)))
        weight, numerator, cost = trial.compile(F(meta['gamma']))
        for offset, polynomial, scale in [(0, weight, cost['metric_scale']), (qsize, numerator, cost['numerator_scale'])]:
            if {flip_mask(m, sites): v for m, v in polynomial.items()} != polynomial:
                raise ValueError('Compiled polynomial is not exactly spin-flip invariant')
            normal = quotient.normal({m: F(v, scale) for m, v in polynomial.items()})
            expected = np.array([float(normal.get(m, 0)) for m in quotient.basis])
            if not np.array_equal(matrix[offset:offset+qsize, j].toarray()[:, 0], expected):
                raise ValueError('Cached metric column differs from exact compilation')
        if (j+1) % 8 == 0:
            print(json.dumps({'phase': 'exact_spinflip_features', 'completed': j+1, 'total': n}), flush=True)
    full_rp = qp + [qsize+i for i in qp] + [2*qsize]
    full_cp = list(range(n)) + [n+i for i in ap] + [n+a+i for i in ap] + [n+2*a]
    if not np.array_equal(rhs, rhs[full_rp]) or (matrix[full_rp, :][:, full_cp] != matrix).nnz:
        raise ValueError('Cached full system violates spin-flip symmetry')
    projection = _average_projection(row_orbits, qsize)
    atom_average = _average_projection(atom_orbits, a).T
    # Identical quotient columns may have different syntactic atom labels.
    reduced_atoms = (projection @ atoms @ atom_average).tocsc()
    reduced_atoms.sum_duplicates(); reduced_atoms.eliminate_zeros(); reduced_atoms.sort_indices()
    seen, keep = {}, []
    for j in range(reduced_atoms.shape[1]):
        sl = slice(reduced_atoms.indptr[j], reduced_atoms.indptr[j+1])
        key = (reduced_atoms.indices[sl].tobytes(), reduced_atoms.data[sl].tobytes())
        if key not in seen:
            seen[key] = j; keep.append(j)
    atom_orbits = [atom_orbits[j] for j in keep]
    atom_average = atom_average[:, keep]
    row_map = block_diag((projection, projection, np.ones((1, 1))), format='csc')
    column_map = block_diag((np.eye(n), atom_average, atom_average, np.ones((1, 1))), format='csc')
    reduced = (row_map @ matrix @ column_map).tocsc()
    reduced.sum_duplicates(); reduced.eliminate_zeros()
    reduced_rhs = row_map @ rhs
    out.mkdir(parents=True, exist_ok=True)
    info = {'source_sha256': input_digest(data), 'gamma': meta['gamma'], 'orbits': orbits,
            'labels': labels, 'atom_orbits': atom_orbits, 'row_orbits': row_orbits,
            'original_shape': list(matrix.shape), 'reduced_shape': list(reduced.shape),
            'nonzeros': reduced.nnz, 'exact_feature_symmetries': n,
            'symmetry_atom_orbits_before_dedup': len(_orbits(ap)),
            'verification_seconds': time.monotonic()-started,
            'scope': 'Generated coefficient proposal reduced by verified spin exchange; orbit weights are averaged. Full rational export required. H6 quotient remains the full fixed-spin function space.'}
    save_npz(out/'model.npz', reduced)
    np.save(out/'rhs.npy', reduced_rhs, allow_pickle=False)
    (out/'prepared.json').write_text(json.dumps(info, indent=2)+'\n')
    print(json.dumps({k: v for k, v in info.items() if k not in ('labels', 'orbits', 'atom_orbits', 'row_orbits')}), flush=True)
    return info


def solve(data, out, time_limit=300, metric_floor=F(9, 10000), numerator_floor=F(1, 10000)):
    import numpy as np
    from scipy.sparse import load_npz
    import highspy
    out = Path(out)
    meta = json.loads((out/'prepared.json').read_text())
    if meta['source_sha256'] != input_digest(data):
        raise ValueError('Reduced model is not bound to bare source')
    matrix = load_npz(out/'model.npz')
    rhs = np.load(out/'rhs.npy', allow_pickle=False)
    n, a = len(meta['orbits']), len(meta['atom_orbits'])
    if matrix.shape != (len(rhs), n+2*a+1):
        raise ValueError('Malformed reduced model')
    rhs[:(len(rhs)-1)//2] *= float(metric_floor)/.001
    lp = highspy.HighsLp()
    lp.num_col_, lp.num_row_ = matrix.shape[1], matrix.shape[0]
    lp.col_cost_ = np.zeros(matrix.shape[1])
    lower = np.zeros(matrix.shape[1]); upper = np.full(matrix.shape[1], highspy.kHighsInf)
    lower[:n] = -highspy.kHighsInf
    lower[-1] = upper[-1] = float(numerator_floor)
    lp.col_lower_, lp.col_upper_ = lower, upper
    lp.row_lower_, lp.row_upper_ = rhs, rhs
    lp.a_matrix_.format_ = highspy.MatrixFormat.kColwise
    lp.a_matrix_.start_, lp.a_matrix_.index_, lp.a_matrix_.value_ = matrix.indptr, matrix.indices, matrix.data
    solver = highspy.Highs()
    for key, value in [('threads', 1), ('solver', 'ipm'), ('presolve', 'on'), ('time_limit', float(time_limit)),
                       ('primal_feasibility_tolerance', 1e-8), ('dual_feasibility_tolerance', 1e-8)]:
        if solver.setOptionValue(key, value) != highspy.HighsStatus.kOk:
            raise ValueError('Native solver rejected option '+key)
    if solver.passModel(lp) == highspy.HighsStatus.kError:
        raise ValueError('Native solver rejected model')
    started = time.monotonic(); solver.run()
    solution = solver.getSolution()
    receipt = {'source_sha256': input_digest(data), 'solver_version': solver.version(),
               'status': solver.modelStatusToString(solver.getModelStatus()), 'solve_seconds': time.monotonic()-started,
               'shape': list(matrix.shape), 'nonzeros': matrix.nnz,
               'metric_floor': str(metric_floor), 'numerator_floor': str(numerator_floor),
               'value_valid': solution.value_valid, 'exact_accepted': False}
    if solution.value_valid:
        x = np.asarray(solution.col_value)
        receipt.update(max_equation_residual=float(max(abs(matrix@x-rhs))), minimum_atom=float(min(x[n:-1])))
        proposal = {'source_sha256': input_digest(data), 'orbits': meta['orbits'],
                    'metric_coefficients': x[:n].tolist(), 'metric_bound': float(metric_floor),
                    'numerator_bound': float(x[-1]), 'gamma': meta['gamma']}
        for name, values in [('weight', x[n:n+a]), ('numerator', x[n+a:n+2*a])]:
            pairs = [(meta['labels'][member], float(value)/len(orbit))
                     for value, orbit in zip(values, meta['atom_orbits']) if abs(value) > 1e-13
                     for member in orbit]
            proposal[name+'_labels'] = [p[0] for p in pairs]
            proposal[name+'_values'] = [p[1] for p in pairs]
        (out/'proposal.json').write_text(json.dumps(proposal, indent=2)+'\n')
        try:
            receipt['exact_receipt'] = export(data, proposal, out/'proof')
            receipt['exact_accepted'] = True
        except ValueError as error:
            receipt['exact_export_rejection'] = str(error)
    (out/'receipt.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps(receipt), flush=True)
    return receipt


if __name__ == '__main__':
    import argparse
    from unittest.mock import patch
    parser = argparse.ArgumentParser()
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--prepared', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--time-limit', type=float, default=300)
    parser.add_argument('--resume', action='store_true')
    args = parser.parse_args()
    with patch('experiments.marginal_determinant_tree.DeterminantOracle.action', side_effect=AssertionError('No states')), \
         patch('experiments.marginal_spin_constructor.spin_states', side_effect=AssertionError('No state generation')), \
         patch('experiments.marginal_polynomial_metric.complete_number_ideals', side_effect=AssertionError('No full-population lift')):
        if not args.resume:
            reduce_cached(args.prepared, args.source, args.out)
        solve(json.loads(args.source.read_text()), args.out, args.time_limit)
