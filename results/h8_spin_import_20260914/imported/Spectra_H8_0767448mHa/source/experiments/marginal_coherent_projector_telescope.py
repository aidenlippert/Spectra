"""Two fixed positive-projector overlap telescopes with four-spin coherences.

A is the particle-hole/spin average of a normalized five-site projector;
Y=(A-RAR*)/2 and T=Y_left-Y_right. Every periodic translated sum vanishes.
A and RAR* are positive contractions, hence ||Y||<=1/2 and ||T||<=1.
"""
from collections import Counter, defaultdict
from fractions import Fraction as F
from experiments.marginal_local_hubbard_block import _exact, _reflection, _sector
from experiments.marginal_charged_projectors import _particlehole, _spinflip

LABELS = {'358,601,-1': {358: 1, 601: -1}, '346,613,1': {346: 1, 613: 1}}


def transformed(vector, action):
    return {action(s)[0]: a*action(s)[1] for s,a in vector.items() if a}


def canonical(vector):
    items=sorted(vector.items()); sign=1 if items[0][1]>0 else -1
    return tuple((s,sign*a) for s,a in items)


def permutation(state, modes):
    occupied = [modes[m] for m in range(len(modes)) if state >> m & 1]
    sign = (-1)**sum(a > b for i, a in enumerate(occupied) for b in occupied[i+1:])
    return sum(1 << m for m in occupied), sign


def spin5(state):
    return permutation(state, [m ^ 1 for m in range(10)])


def hole5(state):
    return 1023 ^ state, (-1)**sum(m+m//2 for m in range(10) if state >> m & 1)


def shift6(state):
    return permutation(state, [(m+2) % 12 for m in range(12)])


def add_outer(matrix, vector, sign):
    for column, b in vector.items():
        for row, a in vector.items():
            matrix[row, column] += sign*a*b


def conjugate(matrix, action):
    result = {}
    for (row, column), value in matrix.items():
        r, a = action(row); c, b = action(column)
        result[r, c] = a*b*value
    return result


def embed(matrix, left):
    result = {}
    for (row, column), value in matrix.items():
        for exterior in range(4):
            r, c = ((row | exterior << 10, column | exterior << 10) if left
                    else (row << 2 | exterior, column << 2 | exterior))
            result[r, c] = value
    return result


def build(vector):
    if (type(vector) is not dict or not 1 <= len(vector) <= 2
            or any(type(s) is not int or not 0 <= s < 1024
                   or type(a) is not int or not 0 < abs(a) <= 10**9 for s, a in vector.items())
            or len({_sector(s, 5) for s in vector}) != 1):
        raise ValueError('One or two bounded integer amplitudes in one five-site spin sector required')
    norm = sum(a*a for a in vector.values())
    images = [vector]
    for action in (hole5, spin5):
        images += [transformed(v, action) for v in list(images)]
    orbit = Counter(canonical(v) for v in images)
    for action in (hole5, spin5):
        actual = Counter()
        for v, multiplicity in orbit.items():
            actual[canonical(transformed(dict(v), action))] += multiplicity
        if actual != orbit:
            raise ValueError('Five-site density orbit is not closed')
    odd = defaultdict(int)
    for v in images:
        if sum(a*a for a in v.values()) != norm:
            raise ValueError('Projector image norm changed')
        add_outer(odd, v, 1)
        add_outer(odd, transformed(v, lambda s: _reflection(s, 5)), -1)
    odd = {key: value for key, value in odd.items() if value}
    denominator = 8*norm
    if conjugate(odd, lambda s: _reflection(s, 5)) != {key: -v for key, v in odd.items()}:
        raise ValueError('Five-site source is not reflection odd')
    left, right = embed(odd, True), embed(odd, False)
    if conjugate(left, shift6) != right:
        raise ValueError('Fermionic translation does not give the right embedding')
    telescope = {key: left.get(key, 0)-right.get(key, 0) for key in left.keys() | right.keys()}
    telescope = {key: value for key, value in telescope.items() if value}
    for (row, column), value in telescope.items():
        if telescope.get((column, row)) != value or _sector(row, 6) != _sector(column, 6):
            raise ValueError('Telescope is not Hermitian or changes spin numbers')
    for action in (_particlehole, _spinflip, lambda s: _reflection(s, 6)):
        if conjugate(telescope, action) != telescope:
            raise ValueError('Six-site telescope breaks a required symmetry')
    # Sum all six fermionic cyclic translates; verify cancellation on the full Fock space.
    summed = defaultdict(int); current = telescope
    for _ in range(6):
        for key, value in current.items():
            summed[key] += value
        current = conjugate(current, shift6)
    if current != telescope or any(summed.values()):
        raise ValueError('Periodic telescope identity failed')
    return denominator, odd, telescope


def coefficients(source):
    if type(source) is not dict or not 1 <= len(source) <= 2:
        raise ValueError('One or two coherent-projector components required')
    if any(type(key) is not str or key not in LABELS for key in source):
        raise ValueError('Canonical coherent-projector label required')
    terms={key: _exact(value) for key,value in source.items()}
    terms={key: value for key,value in terms.items() if value}
    if not terms:
        raise ValueError('Nonzero exact coherent-projector telescope required')
    return terms


def actions(terms):
    result=[{} for _ in range(4096)]
    for label,value in terms.items():
        denominator,_,matrix=build(LABELS[label])
        for (row,column),numerator in matrix.items():
            image=result[column]
            image[row]=image.get(row,F(0))+value*F(numerator,denominator)
    return [{s:a for s,a in image.items() if a} for image in result]
