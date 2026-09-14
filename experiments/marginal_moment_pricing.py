"""Coefficient-moment pricing without a materialized positivity atom matrix."""
from fractions import Fraction as F
from itertools import combinations
from math import comb
import heapq
import json
from pathlib import Path
import time

from experiments.marginal_joint_coefficient_constructor import prepare_features, input_digest, export
from experiments.marginal_joint_spinflip import flip_mask, flip_label, _orbits, _average_projection


class MomentDictionary:
    def __init__(self, data, gamma, feature_degree=4, proof_degree=6):
        import numpy as np
        from scipy.sparse import csc_matrix, vstack
        self.ring, self.quotient, self.orbits, weight, numerator, means = prepare_features(data, F(gamma), feature_degree=feature_degree, quotient_degree=proof_degree)
        ring, quotient = self.ring, self.quotient
        self.sites, self.modes = ring.sites, ring.modes
        self.proof_degree, self.localizer_degree = proof_degree, proof_degree-2
        qp = [quotient.index[flip_mask(m, ring.sites)] for m in quotient.basis]
        self.projection = _average_projection(_orbits(qp), len(qp))
        self.qrows = self.projection.shape[0]
        for polynomial in weight+numerator:
            if {flip_mask(m, ring.sites): v for m, v in polynomial.items()} != polynomial:
                raise ValueError('Exact coefficient target is not spin-exchange invariant')
        def columns(polynomials):
            ri, ci, values = [], [], []
            count = 0
            for j, polynomial in enumerate(polynomials):
                count = j+1
                for mask, value in polynomial.items():
                    ri.append(quotient.index[mask]); ci.append(j); values.append(float(value))
                if len(values) > 12000000:
                    raise ValueError('Coefficient moment-map nonzero budget exceeded')
            return self.projection @ csc_matrix((values, (ri, ci)), shape=(len(quotient.basis), count))
        self.metric = vstack((columns(weight), columns(numerator), csc_matrix([[float(v) for v in means]])), format='csc')
        self.one = (self.projection @ np.array([int(m == 0) for m in quotient.basis])).ravel()
        self.monomials = [sum(1 << i for i in bits) for d in range(proof_degree+1) for bits in combinations(range(ring.modes), d)
                          if all(sum(bool((1 << i)&spin) for i in bits) <= ring.target for spin in ring.spin_masks)]
        self.index = {m: i for i, m in enumerate(self.monomials)}
        # Zero is an extra numerical moment slot for overpopulation monomials.
        self.zero = len(self.monomials)
        polynomials = (quotient.normal({m: F(1)}) for m in self.monomials)
        self.moment_map = columns(polynomials).tocsc()
        self.localizer = {0: -1, **{3 << (2*i): 1 for i in range(ring.sites)}}
        self.charges = [(j, [self.index.get(m | (3 << (2*i)), self.zero) for i in range(ring.sites)])
                        for j, m in enumerate(self.monomials) if m.bit_count() <= self.localizer_degree]
        self.supports = []
        for degree in range(proof_degree+1):
            for bits in combinations(range(ring.modes), degree):
                required = sum(1 << i for i in bits)
                subsets = [0]
                for bit in bits:
                    subsets += [m | (1 << bit) for m in subsets]
                admitted = []
                for j, occupied in enumerate(subsets):
                    count = 1
                    for spin in ring.spin_masks:
                        available = ring.sites-(required&spin).bit_count()
                        needed = ring.target-(occupied&spin).bit_count()
                        count *= comb(available, needed) if 0 <= needed <= available else 0
                    if not count:
                        continue
                    # Canonical spin-exchange atom; no fixed-number dedup needed.
                    if (required, occupied) > (flip_mask(required, ring.sites), flip_mask(occupied, ring.sites)):
                        continue
                    if count > 1:
                        admitted.append((j, ('positive', required, occupied)))
                    if degree <= self.localizer_degree:
                        admitted.append((j, ('charge', required, occupied)))
                if admitted:
                    self.supports.append((degree, [self.index.get(m, self.zero) for m in subsets], admitted))
        self.stats = {'coefficient_rows': self.qrows, 'joint_quotient_dimension': len(quotient.basis),
                      'total_degree': quotient.total_degree, 'charge_feature_degree': feature_degree, 'moment_monomials': len(self.monomials),
                      'moment_map_nonzeros': self.moment_map.nnz,
                      'candidate_labels': sum(len(a) for _, _, a in self.supports),
                      'materialized_all_atom_matrix': False}

    def labels(self, degree=6):
        return [label for d, _, admitted in self.supports if d <= degree for _, label in admitted]

    def column(self, label):
        import numpy as np
        family, required, occupied = label
        polynomial = {occupied: F(1)}
        for i in range(self.modes):
            if (required ^ occupied) & (1 << i):
                polynomial = self.ring.multiply(polynomial, {0: 1, 1 << i: -1})
        if family == 'charge':
            polynomial = self.ring.multiply(polynomial, self.localizer)
        elif family != 'positive':
            raise ValueError('Unsupported atom family')
        normal = self.quotient.normal(polynomial)
        # P averages the coefficient involution, thus represents an orbit average.
        return self.projection @ np.array([float(normal.get(m, 0)) for m in self.quotient.basis])

    def price(self, dual, active, batch):
        import numpy as np
        # Each block of dual acts on P normal(n_M). Pricing atom columns -A
        # therefore negates the finite-difference functional value.
        dual = np.asarray(dual)
        values = self.moment_map.T @ dual[:2*self.qrows].reshape(2, self.qrows).T
        moments = np.vstack((values, np.zeros((1, 2))))
        charge = np.zeros_like(moments)
        for j, targets in self.charges:
            charge[j] = -moments[j]+moments[targets].sum(axis=0)
        heaps = [[], []]; checked = 0; maximum = 0.
        for degree, indices, admitted in self.supports:
            tables = [moments[indices].T.copy()]
            if degree <= self.localizer_degree:
                tables.append(charge[indices].T.copy())
            for table in tables:
                for bit in range(degree):
                    view = table.reshape(2, -1, 2, 1 << bit)
                    view[:, :, 0, :] -= view[:, :, 1, :]
            for index, label in admitted:
                table = tables[0 if label[0] == 'positive' else 1]
                for block in (0, 1):
                    if (block, label) in active:
                        continue
                    score = -float(table[block, index]); checked += 1
                    maximum = max(maximum, score)
                    if score > 1e-8:
                        item = (score, label)
                        if len(heaps[block]) < batch:
                            heapq.heappush(heaps[block], item)
                        elif item > heaps[block][0]:
                            heapq.heapreplace(heaps[block], item)
        return [(b, label) for b in (0, 1) for _, label in sorted(heaps[b], reverse=True)], maximum, checked


def construct(data, gamma, out, batch=32, max_rounds=80, time_limit=180, feature_degree=4, proof_degree=6):
    import numpy as np
    from scipy.sparse import csc_matrix, hstack, eye
    from scipy.optimize import linprog
    started = time.monotonic(); dictionary = MomentDictionary(data, gamma, feature_degree=feature_degree, proof_degree=proof_degree)
    prepared_seconds = time.monotonic()-started
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    rows, n = dictionary.metric.shape
    rhs = np.r_[.0009*dictionary.one, .0001*dictionary.one, 1.]
    active = [(b, label) for b in (0, 1) for label in dictionary.labels(2)]
    active_set = set(active); matrix = csc_matrix((rows, 0))
    built = 0; identity = eye(rows, format='csc'); history = []; accepted = False
    search_started = time.monotonic(); reason = 'Round budget reached'
    for iteration in range(max_rounds):
        remaining = time_limit-(time.monotonic()-search_started)
        if remaining <= 0:
            reason = 'Search time budget reached'; break
        fresh = []
        for block, label in active[built:]:
            column = np.zeros(rows)
            column[block*dictionary.qrows:(block+1)*dictionary.qrows] = -dictionary.column(label)
            fresh.append(column)
        if fresh:
            matrix = hstack((matrix, csc_matrix(np.array(fresh).T)), format='csc')
        built = len(active)
        master = hstack((dictionary.metric, matrix, identity, -identity), format='csc')
        result = linprog(np.r_[np.zeros(n+len(active)), np.ones(2*rows)], A_eq=master, b_eq=rhs,
                         bounds=[(None, None)]*n+[(0, None)]*(len(active)+2*rows), method='highs-ds',
                         options={'time_limit': remaining, 'primal_feasibility_tolerance': 1e-8, 'dual_feasibility_tolerance': 1e-8})
        record = {'round': iteration, 'active_atoms': len(active), 'master_columns': master.shape[1],
                  'master_nonzeros': master.nnz, 'status': result.message,
                  'search_seconds': time.monotonic()-search_started}
        if not result.success:
            reason = 'Restricted solver did not finish'; history.append(record); break
        record['phase_one_l1'] = float(result.fun)
        if result.fun < 1e-7:
            proposal = {'source_sha256': input_digest(data), 'gamma': str(F(gamma)), 'orbits': dictionary.orbits,
                        'metric_coefficients': result.x[:n].tolist(), 'metric_bound': .0009, 'numerator_bound': .0001}
            for block, name in enumerate(('weight', 'numerator')):
                pairs = []
                for (b, label), value in zip(active, result.x[n:n+len(active)]):
                    if b != block or abs(value) <= 1e-13:
                        continue
                    orbit = sorted({label, tuple(flip_label(label, dictionary.sites))})
                    pairs.extend((list(member), float(value)/len(orbit)) for member in orbit)
                proposal[name+'_labels'] = [p[0] for p in pairs]
                proposal[name+'_values'] = [p[1] for p in pairs]
            (out/'proposal.json').write_text(json.dumps(proposal, indent=2)+'\n')
            try:
                record['exact_receipt'] = export(data, proposal, out/'proof')
                accepted = True; reason = 'Exact certificate accepted'
            except ValueError as error:
                reason = 'Exact export rejected'; record['rejection'] = str(error)
            history.append(record); print(json.dumps(record), flush=True); break
        fresh, maximum, checked = dictionary.price(result.eqlin.marginals, active_set, batch)
        record.update(maximum_dual_violation=maximum, candidates_priced=checked, added_atoms=len(fresh))
        history.append(record); print(json.dumps(record), flush=True)
        (out/'history.json').write_text(json.dumps(history, indent=2)+'\n')
        if not fresh:
            reason = 'No numerically improving column; no infeasibility proof'; break
        active.extend(fresh); active_set.update(fresh)
    receipt = dict(dictionary.stats, exact_accepted=accepted, reason=reason, source_sha256=input_digest(data),
                   preparation_seconds=prepared_seconds, elapsed_seconds_including_export=time.monotonic()-started,
                   active_atom_columns=len(active), rounds=len(history),
                   scope='Bare-H generated metric and moment finite-difference pricing; only selected atom columns are materialized. Candidate support labels and bounded-degree moment map remain. No inherited certificates, atom directions, physical configurations, or full atom matrix. Exact full export is acceptance.')
    (out/'receipt.json').write_text(json.dumps(receipt, indent=2)+'\n')
    (out/'history.json').write_text(json.dumps(history, indent=2)+'\n')
    print(json.dumps(receipt), flush=True)
    return receipt


if __name__ == '__main__':
    import argparse
    from unittest.mock import patch
    parser = argparse.ArgumentParser(); parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--gamma', required=True); parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--time-limit', type=float, default=180)
    parser.add_argument('--feature-degree', type=int, default=4)
    parser.add_argument('--proof-degree', type=int, default=6); args = parser.parse_args()
    with patch('experiments.marginal_determinant_tree.DeterminantOracle.action', side_effect=AssertionError('No states')), \
         patch('experiments.marginal_spin_constructor.spin_states', side_effect=AssertionError('No state generation')), \
         patch('experiments.marginal_polynomial_metric.complete_number_ideals', side_effect=AssertionError('No full-population lift')), \
         patch('experiments.marginal_joint_coefficient_constructor.prepare', side_effect=AssertionError('No full atom dictionary')):
        construct(json.loads(args.source.read_text()), args.gamma, args.out, time_limit=args.time_limit, feature_degree=args.feature_degree, proof_degree=args.proof_degree)
