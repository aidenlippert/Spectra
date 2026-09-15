"""Exact counterexample: a fixed-N 2-RDM does not determine all observables or entropy.

Run with standard-library Python. Determinants 7 and 56 occupy disjoint triples
of six fermionic modes. Their equal positive/negative superpositions and their
incoherent mixture have identical one- and two-body reduced density matrices.
"""
from fractions import Fraction as F
from itertools import product
import json


def apply(state, word):
    phase = 1
    for create, mode in reversed(word):
        if ((state >> mode) & 1) == create:
            return None, 0
        phase *= (-1) ** ((state & ((1 << mode) - 1)).bit_count())
        state ^= 1 << mode
    return state, phase


def expect(vector, word):
    norm = sum(a*a for a in vector.values())
    return sum((F(a*phase*vector.get(target, 0), norm)
                for source, a in vector.items()
                for target, phase in [apply(source, word)] if phase), F(0))


def main():
    plus, minus = {7: 1, 56: 1}, {7: 1, 56: -1}
    assert all(s.bit_count() == 3 for s in plus)
    checked = 0
    for degree in (1, 2):
        for indices in product(range(6), repeat=2*degree):
            word = tuple([(1, i) for i in indices[:degree]] +
                         [(0, i) for i in indices[degree:]])
            mixture = (expect({7: 1}, word) + expect({56: 1}, word))/2
            assert expect(plus, word) == expect(minus, word) == mixture
            checked += 1
    triple = ((1, 0), (1, 1), (1, 2), (0, 5), (0, 4), (0, 3))
    adjoint = tuple((1-create, mode) for create, mode in reversed(triple))
    assert expect(plus, triple) + expect(plus, adjoint) == 1
    assert expect(minus, triple) + expect(minus, adjoint) == -1
    # On the two-determinant support, the pure states have eigenvalues (1,0);
    # the incoherent mixture is diagonal with eigenvalues (1/2,1/2).
    for sign in (1, -1):
        rho = [[F(1,2), F(sign,2)], [F(sign,2), F(1,2)]]
        assert [[sum(rho[i][k]*rho[k][j] for k in range(2))
                 for j in range(2)] for i in range(2)] == rho
        assert rho[0][0] + rho[1][1] == 1
    print(json.dumps({'verified': True, 'fixed_N': 3, 'fermionic_modes': 6,
                      'one_and_two_body_words_checked': checked,
                      'identical_1RDM_and_2RDM': True,
                      'three_body_observable_plus': 1,
                      'three_body_observable_minus': -1,
                      'pure_entropy': '0', 'mixture_entropy': 'ln(2)',
                      'scope': 'Information loss of the ordinary fixed-N 2-RDM. Does not exclude richer moment hierarchies or reconstruction under additional assumptions.'}, indent=2))


if __name__ == '__main__':
    main()
