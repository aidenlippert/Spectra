"""Joint polynomial metric/positivity discovery from a bare Hamiltonian.

No physical assignments or inherited metric/positive directions are used.
Numerical proposals remain untrusted until the existing exact replay accepts.
Fixed degree, symmetry and coefficient/atom budgets restrict the search.
"""
from fractions import Fraction as F
from itertools import combinations, product
from math import comb
import json
import hashlib
from pathlib import Path
import time

from experiments.marginal_polynomial_metric import JointPolynomial, replay
from experiments.marginal_number_quotient import NumberSliceQuotient, complete_bounded_number_ideals, _add


def input_digest(data):
    return hashlib.sha256(json.dumps(data, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def feature_orbits(sites, degree=4):
    orbits = {}
    for size in range(min(sites, degree)+1):
        for support in combinations(range(sites), size):
            for powers in product((1, 2), repeat=size):
                if sum(powers) > degree:
                    continue
                vector = [0]*sites
                for i, p in zip(support, powers):
                    vector[i] = p
                vector = tuple(vector)
                key = min(vector, vector[::-1])
                orbits.setdefault(key, set()).add(vector)
    return [sorted(orbits[key]) for key in sorted(orbits)]


def vanishing_feature(orbit, sites):
    terms = {}
    for powers in orbit:
        for i in range(sites):
            vector = list(powers)
            vector[i] += 2
            vector = tuple(0 if p == 0 else 1 if p % 2 else 2 for p in vector)
            terms[vector] = terms.get(vector, 0)+1
    return {'denominator': 2, 'terms': [
        {'powers': list(p), 'coefficient': v} for p, v in sorted(terms.items())]}


def sector_mean(polynomial, sites, population):
    def ratio(k):
        if k > population:
            return F(0)
        return F(comb(sites-k, population-k), comb(sites, population))
    spins = [sum(1 << i for i in range(s, 2*sites, 2)) for s in (0, 1)]
    return sum((v*ratio((m&spins[0]).bit_count())*ratio((m&spins[1]).bit_count())
                for m, v in polynomial.items()), F(0))


class CoefficientQuotient:
    def __init__(self, ring, degree, max_rows=4096):
        if type(degree) is not int or degree < 0 or type(max_rows) is not int or max_rows < 1:
            raise ValueError('Nonnegative total degree and positive row budget required')
        self.total_degree = degree
        degree = min(degree, ring.target)
        ranks = [comb(ring.sites, k)-(comb(ring.sites, k-1) if k else 0) for k in range(degree+1)]
        dimension = sum(a*b for i, a in enumerate(ranks) for j, b in enumerate(ranks) if i+j <= self.total_degree)
        if dimension > max_rows:
            raise ValueError('Joint quotient row budget exceeded')
        self.spins = ring.spin_masks
        self.slices = [NumberSliceQuotient(list(range(s, ring.modes, 2)), ring.target, degree) for s in (0, 1)]
        free = [[m for m in q.monomials if m not in q.rows] for q in self.slices]
        self.basis = [a|b for a in free[0] for b in free[1] if (a|b).bit_count() <= self.total_degree]
        if len(self.basis) != dimension:
            raise ValueError('Unexpected graded joint quotient rank')
        self.index = {m: i for i, m in enumerate(self.basis)}
        self.cache = {}

    def normal(self, polynomial):
        output = {}
        for mask, coefficient in polynomial.items():
            if type(mask) is not int or mask < 0 or mask.bit_count() > self.total_degree or mask & ~(self.spins[0] | self.spins[1]):
                raise ValueError('Polynomial exceeds the supported total-degree coefficient space')
            if mask not in self.cache:
                a, _ = self.slices[0].monomial(mask&self.spins[0])
                b, _ = self.slices[1].monomial(mask&self.spins[1])
                self.cache[mask] = {x|y: u*v for x, u in a.items() for y, v in b.items()}
            _add(output, self.cache[mask], coefficient)
        return output


def prepare_features(data, gamma, feature_degree=4, quotient_degree=6):
    if type(quotient_degree) is not int or not 2 <= quotient_degree <= 6:
        raise ValueError('Coefficient total degree must be in2..6')
    if type(feature_degree) is not int or not 0 <= feature_degree <= 4:
        raise ValueError('Generated charge feature degree must be in0..4')
    if set(data) != {'modes', 'particles', 'hamiltonian'}:
        raise ValueError('Bare Hamiltonian/sector input only; no inherited proof data')
    modes, particles = data['modes'], data['particles']
    if type(modes) is not int or not 4 <= modes <= 64 or modes % 4 or type(particles) is not int or particles != modes//2:
        raise ValueError('Even-site half-filled sector within the verifier mode budget required')
    sites = modes//2
    ring = JointPolynomial(dict(data, polynomial_metric=vanishing_feature([[0]*sites], sites)))
    quotient = CoefficientQuotient(ring, quotient_degree)
    orbits = feature_orbits(sites, feature_degree)
    q = [{0: -1, 1 << (2*i): 1, 1 << (2*i+1): 1} for i in range(sites)]
    selected, weights, pivots, means = [], [], {}, []
    for orbit in orbits:
        metric = vanishing_feature(orbit, sites)
        trial = JointPolynomial(dict(data, polynomial_metric=metric))
        polynomial = {m: F(v, 2) for m, v in trial.metric_polynomial(q).items()}
        normal = quotient.normal(polynomial)
        reduced = dict(normal)
        while reduced:
            pivot = min(reduced, key=quotient.index.__getitem__)
            coefficient = reduced[pivot]
            if pivot not in pivots:
                pivots[pivot] = {m: v/coefficient for m, v in reduced.items()}
                selected.append(orbit)
                weights.append(normal)
                dimension = comb(sites, ring.target)
                # Every feature vanishes on P, hence E_Q=dimension/(dimension-1)*E_sector.
                means.append(F(dimension, dimension-1)*sector_mean(polynomial, sites, ring.target))
                break
            _add(reduced, pivots[pivot], -coefficient)
    print(json.dumps({'phase': 'feature_basis', 'generated_orbits': len(orbits),
                      'independent_features': len(selected)}), flush=True)
    numerators = []
    for i, orbit in enumerate(selected):
        trial = JointPolynomial(dict(data, polynomial_metric=vanishing_feature(orbit, sites)))
        _, polynomial, cost = trial.compile(gamma)
        numerators.append(quotient.normal({m: F(v, cost['numerator_scale']) for m, v in polynomial.items()}))
        if i % 8 == 7:
            print(json.dumps({'phase': 'feature_compilation', 'completed': i+1, 'total': len(selected)}), flush=True)
    return ring, quotient, selected, weights, numerators, means


def prepare(data, gamma, max_columns=100000, max_nonzeros=12000000):
    from scipy.sparse import csc_matrix
    ring, quotient, selected, weights, numerators, means = prepare_features(data, gamma)
    sites = ring.sites
    labels, ri, ci, values = [], [], [], []
    localizer = {0: -1, **{3 << (2*i): 1 for i in range(sites)}}
    for degree in range(7):
        for bits in combinations(range(ring.modes), degree):
            required = sum(1 << i for i in bits)
            for assignment in range(1 << degree):
                occupied = sum(1 << i for j, i in enumerate(bits) if assignment & (1 << j))
                count = 1
                for spin in ring.spin_masks:
                    available = sites-(required&spin).bit_count()
                    needed = ring.target-(occupied&spin).bit_count()
                    count *= comb(available, needed) if 0 <= needed <= available else 0
                if not count:
                    continue
                atom = {occupied: F(1)}
                for i in bits:
                    if not occupied & (1 << i):
                        atom = ring.multiply(atom, {0: 1, 1 << i: -1})
                families = ([] if count == 1 else [('positive', atom)])
                if degree <= 4:
                    families.append(('charge', ring.multiply(atom, localizer)))
                for family, polynomial in families:
                    normal = quotient.normal(polynomial)
                    if not normal:
                        continue
                    if len(labels) >= max_columns or len(values)+len(normal) > max_nonzeros:
                        raise ValueError('Positivity dictionary budget exceeded')
                    column = len(labels)
                    labels.append([family, required, occupied])
                    for mask, v in normal.items():
                        ri.append(quotient.index[mask]); ci.append(column); values.append(float(v))
        print(json.dumps({'phase': 'atom_dictionary', 'degree': degree, 'columns': len(labels), 'nonzeros': len(values)}), flush=True)
    matrix = csc_matrix((values, (ri, ci)), shape=(len(quotient.basis), len(labels)))
    return ring, quotient, selected, weights, numerators, means, labels, matrix


def construct(data, gamma, out, time_limit=120, method='highs-ipm'):
    import numpy as np
    from scipy.sparse import csc_matrix, bmat
    from scipy.optimize import linprog
    start = time.monotonic()
    gamma = F(gamma)
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    ring, quotient, orbits, weight, numerator, means, labels, atoms = prepare(data, gamma)
    def columns(polynomials):
        matrix = np.zeros((len(quotient.basis), len(polynomials)))
        for j, polynomial in enumerate(polynomials):
            for m, v in polynomial.items():
                matrix[quotient.index[m], j] = float(v)
        return csc_matrix(matrix)
    one = np.zeros((len(quotient.basis), 1)); one[quotient.index[0], 0] = 1
    matrix = bmat([[columns(weight), -atoms, None, None],
                   [columns(numerator), None, -atoms, csc_matrix(-one)],
                   [csc_matrix([[float(x) for x in means]]), None, None, csc_matrix((1, 1))]], format='csc')
    rhs = np.r_[.001*one[:, 0], np.zeros(len(quotient.basis)), 1.]
    objective = np.zeros(matrix.shape[1]); objective[-1] = -1
    bounds = [(None, None)]*len(orbits)+[(0, None)]*(2*len(labels))+[(None, None)]
    prepared = time.monotonic()
    from scipy.sparse import save_npz
    save_npz(out/'model.npz', matrix)
    np.save(out/'rhs.npy', rhs, allow_pickle=False)
    metadata = {'source_sha256': input_digest(data), 'gamma': str(gamma), 'orbits': orbits,
                'labels': labels, 'quotients': [q.stats for q in quotient.slices],
                'preparation_seconds': prepared-start}
    (out/'prepared.json').write_text(json.dumps(metadata, indent=2)+'\n')
    return solve_prepared(data, gamma, out, time_limit, method)


def solve_prepared(data, gamma, out, time_limit=180, method='highs-ipm', presolve=False):
    import numpy as np
    from scipy.sparse import load_npz
    from scipy.optimize import linprog
    out = Path(out)
    metadata = json.loads((out/'prepared.json').read_text())
    if metadata['source_sha256'] != input_digest(data) or F(metadata['gamma']) != F(gamma):
        raise ValueError('Prepared model is not bound to requested H and threshold')
    matrix = load_npz(out/'model.npz')
    rhs = np.load(out/'rhs.npy', allow_pickle=False)
    orbits, labels = metadata['orbits'], metadata['labels']
    if matrix.shape != (len(rhs), len(orbits)+2*len(labels)+1):
        raise ValueError('Malformed prepared model')
    objective = np.zeros(matrix.shape[1]); objective[-1] = -1
    bounds = [(None, None)]*len(orbits)+[(0, None)]*(2*len(labels))+[(None, None)]
    print(json.dumps({'phase': 'joint_solve', 'rows': matrix.shape[0], 'columns': matrix.shape[1], 'nonzeros': matrix.nnz, 'method': method, 'presolve': presolve}), flush=True)
    started = time.monotonic()
    result = linprog(objective, A_eq=matrix, b_eq=rhs, bounds=bounds, method=method,
                     options={'time_limit': time_limit, 'presolve': presolve, 'disp': True})
    receipt = {'success': bool(result.success), 'status': result.message, 'target_lower': str(F(gamma)),
               'source_sha256': input_digest(data), 'preparation_seconds': metadata['preparation_seconds'],
               'solve_seconds': time.monotonic()-started, 'iterations': result.nit,
               'method': method, 'presolve': presolve,
               'rows': matrix.shape[0], 'columns': matrix.shape[1], 'nonzeros': matrix.nnz,
               'metric_features': len(orbits), 'atom_columns': len(labels),
               'quotients': metadata['quotients'], 'configuration_evaluations': 0,
               'scope': 'Bare-Hamiltonian numerical joint discovery; generated reflection-tied D-times-degree4 charge basis and complete nonsingleton degree6/charge-localizer degree4 dictionary. No inherited coefficients, selected atoms, physical assignments or actions. Fixed degree and resource budgets restrict the search; exact export required.'}
    if result.success:
        n = len(orbits); a = len(labels)
        receipt.update(numerator_bound=float(result.x[-1]), max_equation_residual=float(max(abs(matrix@result.x-rhs))),
                       minimum_atom_coefficient=float(min(result.x[n:-1])))
        proposal = {'source_sha256': input_digest(data), 'orbits': orbits, 'metric_coefficients': result.x[:n].tolist(),
                    'metric_bound': .001, 'numerator_bound': float(result.x[-1]), 'gamma': str(F(gamma))}
        for name, values in [('weight', result.x[n:n+a]), ('numerator', result.x[n+a:n+2*a])]:
            pairs = [(label, float(value)) for label, value in zip(labels, values) if abs(value) > 1e-13]
            proposal[name+'_labels'] = [p[0] for p in pairs]
            proposal[name+'_values'] = [p[1] for p in pairs]
        (out/'proposal.json').write_text(json.dumps(proposal, indent=2)+'\n')
    (out/'numerical_receipt.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps(receipt), flush=True)
    return receipt


def export(data, proposal, out, metric_precision=10**9, proof_denominator=10**12):
    if set(data) != {'modes', 'particles', 'hamiltonian'}:
        raise ValueError('Bare Hamiltonian/sector input only')
    if proposal.get('source_sha256') != input_digest(data):
        raise ValueError('Proposal is not bound to this bare Hamiltonian')
    if type(metric_precision) is not int or not 1<=metric_precision<=5*10**13:
        raise ValueError('Metric precision must preserve the verifier denominator bound')
    if type(proof_denominator) is not int or not 1<=proof_denominator<=10**16:
        raise ValueError('Proof denominator must preserve the verifier bound')
    metric_denominator, denominator = 2*metric_precision, proof_denominator
    terms = {}
    sites = data['modes']//2
    for orbit, coefficient in zip(proposal['orbits'], proposal['metric_coefficients']):
        integer = round(coefficient*metric_precision)
        for term in vanishing_feature(orbit, sites)['terms']:
            key = tuple(term['powers'])
            terms[key] = terms.get(key, 0)+integer*term['coefficient']
    certificate = dict(data, kind='joint_polynomial_metric_gap_v1', target_lower=proposal['gamma'],
        polynomial_metric={'denominator': metric_denominator, 'terms': [
            {'powers': list(p), 'coefficient': v} for p, v in sorted(terms.items()) if v]})
    ring = JointPolynomial(certificate)
    weight, numerator, cost = ring.compile(F(certificate['target_lower']))
    diagnostics = {'scope': 'Exact rational export of bare-H joint proposal, using bounded-degree number identities. No inherited certificate coefficients or directions.', 'parts': {}}
    localizer = {0: -1, **{3 << (2*i): 1 for i in range(sites)}}
    for name, polynomial, scale in [('weight', weight, cost['metric_scale']), ('numerator', numerator, cost['numerator_scale'])]:
        bound = round(proposal['metric_bound' if name == 'weight' else 'numerator_bound']*denominator)
        proof = {'denominator': denominator, 'bound': bound,
                 'positive_indicators': [], 'charge_indicators': [], 'number_multipliers': [[], []]}
        represented = {0: F(bound, denominator)}
        negative = []
        for label, value in zip(proposal[name+'_labels'], proposal[name+'_values']):
            if value < 0:
                negative.append(value)
            integer = max(0, round(value*denominator))
            if not integer:
                continue
            family, required, occupied = label
            if family not in ('positive', 'charge'):
                raise ValueError('Unknown positivity atom')
            proof['positive_indicators' if family == 'positive' else 'charge_indicators'].append(
                {'required': required, 'occupied': occupied, 'weight': integer})
            term = {occupied: F(integer, denominator)}
            for i in range(ring.modes):
                if (required^occupied) & (1 << i):
                    term = ring.multiply(term, {0: 1, 1 << i: -1})
            if family == 'charge':
                term = ring.multiply(term, localizer)
            _add(represented, term)
        residual = ring.add({m: F(v, scale) for m, v in polynomial.items()}, {m: -v for m, v in represented.items()})
        remainder, ideals, stats = complete_bounded_number_ideals(ring, residual)
        proof['number_multipliers'] = [[{'mask': m, 'coefficient': round(v*denominator)}
                                       for m, v in sorted(ideal.items()) if round(v*denominator)] for ideal in ideals]
        certificate[name+'_proof'] = proof
        diagnostics['parts'][name] = {'quotients': stats, 'quotient_remainder_l1': str(sum(map(abs, remainder.values()), F(0))),
            'clipped_negative_count': len(negative), 'most_negative_proposal': min(negative, default=0)}
    dimension = comb(sites, ring.target)
    diagnostics['exact_mean_Q_metric'] = str(F(dimension, dimension-1)*sector_mean({m: F(v, cost['metric_scale']) for m, v in weight.items()}, sites, ring.target))
    receipt = replay(certificate)
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    for name, value in [('certificate', certificate), ('receipt', receipt), ('construction', diagnostics)]:
        (out/(name+'.json')).write_text(json.dumps(value, indent=2)+'\n')
    return receipt


if __name__ == '__main__':
    import argparse
    from unittest.mock import patch
    parser = argparse.ArgumentParser()
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--gamma', required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--time-limit', type=float, default=120)
    parser.add_argument('--export-only', action='store_true')
    parser.add_argument('--resume', action='store_true')
    parser.add_argument('--method', choices=['highs-ipm', 'highs-ds'], default='highs-ipm')
    args = parser.parse_args()
    with patch('experiments.marginal_determinant_tree.DeterminantOracle.action', side_effect=AssertionError('No physical states')), \
         patch('experiments.marginal_spin_constructor.spin_states', side_effect=AssertionError('No sector generation')), \
         patch('experiments.marginal_polynomial_metric.complete_number_ideals', side_effect=AssertionError('No full-population lift')):
        data = json.loads(args.source.read_text())
        success = True
        if not args.export_only:
            if args.resume:
                success = solve_prepared(data, args.gamma, args.out, args.time_limit, args.method)['success']
            else:
                success = construct(data, args.gamma, args.out, args.time_limit, args.method)['success']
        proposal_path = args.out/'proposal.json'
        if success and proposal_path.exists():
            proposal = json.loads(proposal_path.read_text())
            if F(proposal['gamma']) != F(args.gamma):
                raise ValueError('Proposal threshold differs from requested target')
            if proposal['numerator_bound'] > 0:
                print(json.dumps(export(data, proposal, args.out/'proof')), flush=True)
