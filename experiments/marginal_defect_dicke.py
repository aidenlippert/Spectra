"""Exact local-defect/Dicke algebra without configuration enumeration.

An atom (occupation_tuple, k) fixes named pairs to local occupations 0,1,2,3
(empty,L,R,LR) and sums all singly occupied spectator pairs with k rights.
The Fock gauge orders modes L0,R0,L1,R1,...; canonical CAR words are mapped
into that gauge before applying their Jordan-Wigner signs. Atoms overlap
and are not a linearly independent basis. Zero norm, not dictionary zero,
is the valid test for equality of represented physical vectors.
"""
from fractions import Fraction as F
from functools import lru_cache
from math import comb, lcm
import json
from pathlib import Path
import time

from experiments.marginal_collective import hopping_polynomial
from experiments.marginal_general_schur import norm_bound, complement_lower, pivots, upper_rayleigh
from experiments.marginal_symbolic import add, canonical, decode, hermitian, scale


def combine(*pieces):
    result = {}
    for weight, state in pieces:
        if not weight:
            continue
        for atom, value in state.items():
            result[atom] = result.get(atom, F(0)) + weight * value
    return {atom: value for atom, value in result.items() if value}


class DefectDicke:
    def __init__(self, pairs):
        if type(pairs) is not int or pairs < 2:
            raise ValueError('At least two pairs required')
        self.pairs = pairs
        # Cache belongs to this instance; do not retain past model instances.
        self.atom_overlap = lru_cache(maxsize=200000)(self._atom_overlap)
        self.coordinate_sectors = {}
        self.atom_coordinate_cache = {}

    def atom(self, occupations=(), rights=0):
        occupations = tuple(sorted(occupations))
        if (type(rights) is not int or len({i for i, _ in occupations}) != len(occupations)
                or any(type(i) is not int or not 0 <= i < self.pairs or type(o) is not int or not 0 <= o < 4
                       for i, o in occupations)):
            raise ValueError('Invalid defect atom')
        if not 0 <= rights <= self.pairs - len(occupations):
            return {}
        return {(occupations, rights): F(1)}

    def reference(self):
        return [self.atom(rights=k) for k in range(self.pairs + 1)]

    def _atom_overlap(self, left, right):
        a, ka = left
        b, kb = right
        a, b = dict(a), dict(b)
        for i in a.keys() & b.keys():
            if a[i] != b[i]:
                return 0
        for i in a.keys() - b.keys():
            if a[i] not in (1, 2):
                return 0
            kb -= a[i] == 2
        for i in b.keys() - a.keys():
            if b[i] not in (1, 2):
                return 0
            ka -= b[i] == 2
        n = self.pairs - len(a.keys() | b.keys())
        return comb(n, ka) if ka == kb and 0 <= ka <= n else 0

    def atom_coordinates(self, atom):
        """Exact kernel Gram-LDL coordinates of an encountered atom.

        New orthogonal directions have zero coordinates for all older atoms,
        so cached short coordinate vectors remain valid as a sector grows.
        """
        if atom in self.atom_coordinate_cache:
            return self.atom_coordinate_cache[atom]
        if len(self.atom_coordinate_cache) >= 200000:
            raise ValueError('Implicit atom-coordinate budget exhausted')
        charge = tuple((i, o) for i, o in atom[0] if o in (0, 3))
        rights = atom[1] + sum(bool(o & 2) for _, o in atom[0])
        key = (charge, rights)
        sector = self.coordinate_sectors.setdefault(key, {'atoms': [], 'lower': [], 'diagonal': []})
        coordinates = []
        diagonal = sector['diagonal']
        for i, basis in enumerate(sector['atoms']):
            overlap = self.atom_overlap(basis, atom) if basis <= atom else self.atom_overlap(atom, basis)
            coordinates.append((overlap - sum(sector['lower'][i][j] * diagonal[j] * coordinates[j]
                                               for j in range(i))) / diagonal[i])
        residual = F(self.atom_overlap(atom, atom)) - sum(d * x * x for d, x in zip(diagonal, coordinates))
        if residual < 0:
            raise ValueError('Negative exact atom-kernel residual')
        if residual:
            sector['atoms'].append(atom)
            sector['lower'].append(coordinates + [F(1)])
            diagonal.append(residual)
            coordinates.append(F(1))
        result = (key, tuple(coordinates))
        self.atom_coordinate_cache[atom] = result
        return result

    def _coordinate_state(self, state):
        denominator = lcm(*(value.denominator for value in state.values()))
        sectors = {}
        for atom, value in state.items():
            key, coordinates = self.atom_coordinates(atom)
            target = sectors.setdefault(key, [])
            target.extend([F(0)] * (len(coordinates) - len(target)))
            coefficient = value.numerator * (denominator // value.denominator)
            for i, x in enumerate(coordinates):
                if x:
                    target[i] += coefficient * x
        return sectors, denominator

    def _coordinate_inner(self, left, right):
        result = F(0)
        for key, coordinates in left[0].items():
            if key in right[0]:
                result += sum(d * a * b for d, a, b in zip(self.coordinate_sectors[key]['diagonal'], coordinates, right[0][key]))
        return result / (left[1] * right[1])

    def inner(self, left, right):
        return self._coordinate_inner(self._coordinate_state(left), self._coordinate_state(right))

    def gram(self, left, right):
        a = [self._coordinate_state(v) for v in left]
        if left is right:
            result = [[F(0)] * len(a) for _ in a]
            for i in range(len(a)):
                for j in range(i + 1):
                    result[i][j] = result[j][i] = self._coordinate_inner(a[i], a[j])
            return result
        b = [self._coordinate_state(v) for v in right]
        return [[self._coordinate_inner(x, y) for y in b] for x in a]

    def apply_word(self, word, state):
        current = state
        for creation, mode in reversed(word):
            if creation not in (0, 1) or type(mode) is not int or not 0 <= mode < 2 * self.pairs:
                raise ValueError('Invalid CAR mode')
            pair = mode % self.pairs
            bit = 1 if mode < self.pairs else 2
            result = {}
            for (occupation_tuple, rights), coefficient in current.items():
                occupations = dict(occupation_tuple)
                choices = [(occupations[pair], rights)] if pair in occupations else [(1, rights), (2, rights - 1)]
                for occupancy, remaining in choices:
                    n = self.pairs - len(occupations) - (pair not in occupations)
                    if not 0 <= remaining <= n or bool(occupancy & bit) == bool(creation):
                        continue
                    prefix = pair + sum(o.bit_count() - 1 for i, o in occupations.items() if i < pair)
                    if bit == 2 and occupancy & 1:
                        prefix += 1
                    updated = dict(occupations)
                    updated[pair] = occupancy ^ bit
                    atom = (tuple(sorted(updated.items())), remaining)
                    result[atom] = result.get(atom, F(0)) + coefficient * (-1) ** prefix
            current = {atom: value for atom, value in result.items() if value}
        return current

    def apply_polynomial(self, polynomial, state):
        return combine(*[(coefficient, self.apply_word(word, state)) for word, coefficient in polynomial.items()])

    def reference_action(self, state, hopping=F(1, 5)):
        """H0 action keeps the distinguished-pair set fixed."""
        result = {}
        def put(atom, value):
            if value:
                result[atom] = result.get(atom, F(0)) + value
        for (occupations, k), coefficient in state.items():
            n = self.pairs - len(occupations)
            nl = n - k + sum(bool(o & 1) for _, o in occupations)
            nr = k + sum(bool(o & 2) for _, o in occupations)
            put((occupations, k), coefficient * (comb(nl, 2) + comb(nr, 2)))
            for index, (i, o) in enumerate(occupations):
                if o in (1, 2):
                    updated = occupations[:index] + ((i, 3 - o),) + occupations[index + 1:]
                    put((updated, k), -hopping * coefficient)
            if k < n:
                put((occupations, k + 1), -hopping * coefficient * (k + 1))
            if k:
                put((occupations, k - 1), -hopping * coefficient * (n - k + 1))
        return {atom: value for atom, value in result.items() if value}

    def coupling(self, h):
        h = canonical(h)
        if not hermitian(h):
            raise ValueError('Real Hermitian Hamiltonian required')
        if any(len(word) > 4 or sum(1 if creation else -1 for creation, _ in word) for word in h):
            raise ValueError('Number-conserving degree <=4 required')
        if any(not 0 <= mode < 2 * self.pairs for word in h for _, mode in word):
            raise ValueError('Invalid mode')
        z = self.reference()
        metric = [F(comb(self.pairs, k)) for k in range(self.pairs + 1)]
        h0 = hopping_polynomial(2 * self.pairs, F(1, 5))
        delta = add(h, scale(h0, -1))
        dz = [self.apply_polynomial(delta, v) for v in z]
        projected = self.gram(z, dz)
        w = [combine((1, v), *[(-projected[i][j] / metric[i], z[i]) for i in range(len(z))])
             for j, v in enumerate(dz)]
        if any(x for row in self.gram(z, w) for x in row):
            raise ValueError('Exact coupling orthogonality failed')
        hz = [combine((1, self.reference_action(v)), (1, dv)) for v, dv in zip(z, dz)]
        return {'embedding': z, 'coupling': w, 'metric': metric, 'projected_h': self.gram(z, hz),
                'perturbation_norm_bound': norm_bound(delta, 'hermitian_pairs')}

    def recurrence(self, columns, max_degree=20):
        """Exact scalar Lanczos in the Hilbert space of block columns.

        The zero residual is checked through its exact squared norm; physical
        identities among redundant atoms are therefore retained correctly.
        """
        def inner(a, b):
            return sum(self.inner(x, y) for x, y in zip(a, b))
        q = columns
        previous = [{} for _ in columns]
        polynomial, previous_polynomial = [F(1)], []
        previous_norm = F(1)
        norms = []
        for degree in range(max_degree + 1):
            norm = inner(q, q)
            if norm < 0:
                raise ValueError('Negative exact squared norm')
            if not norm:
                return {'annihilator': polynomial, 'lanczos_norms': norms,
                        'zero_residual_norm': '0'}
            if degree == max_degree:
                raise ValueError('Exact recurrence degree budget exhausted')
            hq = [self.reference_action(v) for v in q]
            alpha = inner(q, hq) / norm
            beta = norm / previous_norm if degree else F(0)
            next_q = [combine((1, hv), (-alpha, v), (-beta, pv)) for hv, v, pv in zip(hq, q, previous)]
            next_polynomial = [F(0)] + polynomial
            for i, value in enumerate(polynomial):
                next_polynomial[i] -= alpha * value
            for i, value in enumerate(previous_polynomial):
                next_polynomial[i] -= beta * value
            norms.append(str(norm))
            previous, q = q, next_q
            previous_polynomial, polynomial = polynomial, next_polynomial
            previous_norm = norm
        raise AssertionError('Unreachable')

    def closure(self, columns, max_dimension=64):
        """Select exact independent H0-Krylov vectors by Gram LDL.

        Raw atoms are redundant, so coefficient-row rank would be incorrect.
        Positive Schur residual norms admit vectors; zero norms reject them.
        Every admitted vector's H0 image enters the finite work queue.
        """
        basis, triangular, diagonal = [], [], []
        queue = list(columns)
        cursor = 0
        while cursor < len(queue):
            v = queue[cursor]
            cursor += 1
            overlaps = [self.inner(b, v) for b in basis]
            row = []
            for i, overlap in enumerate(overlaps):
                row.append((overlap - sum(triangular[i][j] * diagonal[j] * row[j] for j in range(i))) / diagonal[i])
            residual = self.inner(v, v) - sum(d * value * value for d, value in zip(diagonal, row))
            if residual < 0:
                raise ValueError('Negative exact Gram residual')
            if not residual:
                continue
            if len(basis) >= max_dimension:
                raise ValueError('Implicit closure dimension budget exhausted')
            basis.append(v)
            triangular.append(row + [F(1)])
            diagonal.append(residual)
            queue.append(self.reference_action(v))
        return basis


def prepare(h, resolvent=True):
    model = DefectDicke(5)
    data = model.coupling(h)
    w = data['coupling']
    data['leakage'] = model.gram(w, w)
    data['norm_method'] = 'hermitian_pairs'
    data['complement_lower'] = complement_lower()
    if resolvent:
        data.update(model.recurrence(w, max_degree=12))
        moments = []
        power = w
        for _ in range(len(data['annihilator']) - 1):
            moments.append(model.gram(w, power))
            power = [model.reference_action(v) for v in power]
        data['moments'] = moments
    data['implicit_atom_count'] = sum(len(v) for v in w)
    data['max_active_pairs'] = max((len(atom[0]) for v in w for atom in v), default=0)
    return data


def replay(certificate):
    if certificate.get('kind') != 'general_perturbation_schur_v1' or certificate.get('modes') != 10 or certificate.get('particles') != 5:
        raise ValueError('Expected ten-mode general Schur certificate')
    if certificate.get('norm_method') != 'hermitian_pairs' or type(certificate.get('resolvent')) is not bool:
        raise ValueError('Expected paired-norm certificate with explicit resolvent flag')
    h = decode(certificate['hamiltonian'], 10, 4)
    data = prepare(h, certificate['resolvent'])
    lower = F(certificate['lower'])
    checked = pivots(data, lower)
    if checked is None:
        raise ValueError('Nonpositive exact Schur certificate')
    upper = upper_rayleigh(h, certificate['independent_upper'])
    if lower > upper:
        raise ValueError('Inconsistent certified interval')
    return {'lower': str(lower), 'upper': str(upper), 'width': str(upper - lower),
            'width_float': float(upper - lower), 'schur_pivots': [str(x) for x in checked],
            'implicit_atom_count': data['implicit_atom_count'], 'max_active_pairs': data['max_active_pairs'],
            'annihilator': [str(x) for x in data.get('annihilator', [])],
            'scope': 'Lower Schur construction uses implicit defect-Dicke algebra; existing upper witness remains explicit'}


def reference_complement(pairs):
    """All complementary spins lie above j=pairs/2-1 (see size-transfer proof)."""
    from experiments.marginal_sector_reference import bracket
    if type(pairs) is not int or pairs < 2:
        raise ValueError('At least two pairs required')
    return F(bracket(2 * pairs, F(1, 5), 1)['lower'])


def enlarged_workspace(h, pairs=5):
    model = DefectDicke(pairs)
    if pairs + 1 > 64:
        raise ValueError('Reference dimension exceeds retained-space budget')
    base = model.coupling(h)
    base['complement_lower'] = reference_complement(pairs)
    r = model.closure(base['coupling'], max_dimension=64 - len(base['embedding']))
    u = base['embedding'] + r
    if any(x for row in model.gram(base['embedding'], r) for x in row):
        raise ValueError('Implicit closure is not orthogonal to reference space')
    delta = add(h, scale(hopping_polynomial(2 * pairs, F(1, 5)), -1))
    hu = [combine((1, model.reference_action(v)), (1, model.apply_polynomial(delta, v))) for v in u]
    return {'model': model, 'basis': u, 'action': hu, 'delta': delta, 'base': base}


def workspace_data(workspace):
    from experiments.marginal_enlarged_schur import solve_positive
    if 'data' in workspace:
        return workspace['data']
    model, u, hu = (workspace[k] for k in ('model', 'basis', 'action'))
    g = model.gram(u, u)
    a = model.gram(u, hu)
    coefficients = solve_positive(g, a)
    # L^T L = (HU)^T HU - A^T G^-1 A, without expanding projected vectors.
    second = model.gram(hu, hu)
    n = len(u)
    leakage = [[second[i][j] - sum(a[k][i] * coefficients[k][j] for k in range(n))
                for j in range(n)] for i in range(n)]
    if any(leakage[i][j] for i in range(len(workspace['base']['embedding'])) for j in range(n)):
        raise ValueError('Initial coupling was not retained exactly')
    if any(a[i][j] != a[j][i] or leakage[i][j] != leakage[j][i] for i in range(n) for j in range(n)):
        raise ValueError('Nonsymmetric exact projected data')
    data = {'metric': g, 'projected_h': a, 'leakage': leakage,
            'norm_bound': workspace['base']['perturbation_norm_bound'], 'complement_lower': workspace['base']['complement_lower'],
            'retained_dimension': n, 'basis_atom_count': sum(len(v) for v in u),
            'action_atom_count': sum(len(v) for v in hu),
            'max_active_pairs': max((len(atom[0]) for v in hu for atom in v), default=0)}
    workspace['data'] = data
    return data


def enlarged_prepare(h):
    return workspace_data(enlarged_workspace(h))


def targeted_workspace(workspace, coefficients):
    from experiments.marginal_enlarged_schur import solve_positive
    model, u, hu = (workspace[k] for k in ('model', 'basis', 'action'))
    if type(coefficients) is not list or len(coefficients) != len(u) or any(type(x) is not int for x in coefficients) or not any(coefficients):
        raise ValueError('Expected nonzero integer targeting coordinates')
    data = workspace_data(workspace)
    rhs = [[sum(row[j] * coefficients[j] for j in range(len(u)))] for row in data['projected_h']]
    projection = solve_positive(data['metric'], rhs)
    seed = combine(*[(c, v) for c, v in zip(coefficients, hu)], *[(-row[0], v) for row, v in zip(projection, u)])
    if any(model.inner(v, seed) for v in u):
        raise ValueError('Targeted leakage projection failed')
    r = model.closure([seed], max_dimension=64 - len(u))
    if any(x for row in model.gram(u, r) for x in row):
        raise ValueError('Targeted closure is not orthogonal')
    hr = [combine((1, model.reference_action(v)), (1, model.apply_polynomial(workspace['delta'], v))) for v in r]
    return {'model': model, 'basis': u + r, 'action': hu + hr,
            'delta': workspace['delta'], 'base': workspace['base'], 'previous_dimension': len(u)}


def enlarged_replay(certificate):
    from experiments.marginal_enlarged_schur import pivots as enlarged_pivots
    if (certificate.get('kind') != 'enlarged_coupling_schur_v1' or certificate.get('modes') != 10
            or certificate.get('particles') != 5 or certificate.get('rounds', 1) != 1
            or certificate.get('resolvent', False) is not False
            or certificate.get('target_coefficients') is not None):
        raise ValueError('Expected one-round scalar enlarged certificate')
    h = decode(certificate['hamiltonian'], 10, 4)
    data = enlarged_prepare(h)
    lower = F(certificate['lower'])
    checked = enlarged_pivots(data, lower)
    if checked is None:
        raise ValueError('Nonpositive exact enlarged Schur certificate')
    upper = upper_rayleigh(h, certificate['independent_upper'])
    if lower > upper:
        raise ValueError('Inconsistent interval')
    return {'lower': str(lower), 'upper': str(upper), 'width': str(upper - lower),
            'width_float': float(upper - lower), 'retained_dimension': data['retained_dimension'],
            'schur_pivots': [str(x) for x in checked], 'basis_atom_count': data['basis_atom_count'],
            'action_atom_count': data['action_atom_count'], 'max_active_pairs': data['max_active_pairs'],
            'scope': 'One-round enlarged lower proof uses implicit defect-Dicke algebra; upper witness remains explicit'}


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify', required=True)
    args = parser.parse_args()
    started = time.monotonic()
    certificate = json.loads(Path(args.verify).read_text())
    verifier = enlarged_replay if certificate.get('kind') == 'enlarged_coupling_schur_v1' else replay
    result = verifier(certificate)
    result['elapsed_seconds'] = time.monotonic() - started
    print(json.dumps(result, indent=2))
