import unittest
import numpy as np
from fractions import Fraction

from .model import fixture, mpo, numerical_mpo, trace_seed, seed_squared_norm


def car_word_matrix(word, modes):
    n = 1 << modes; out = np.zeros((n, n))
    for col in range(n):
        state = col; amp = 1
        # product is written left-to-right, so act rightmost first
        for create, mode in reversed(word):
            if create:
                if state & (1 << mode): amp = 0; break
                amp *= -1 if (state & ((1 << mode) - 1)).bit_count() % 2 else 1
                state |= 1 << mode
            else:
                if not state & (1 << mode): amp = 0; break
                amp *= -1 if (state & ((1 << mode) - 1)).bit_count() % 2 else 1
                state ^= 1 << mode
        if amp: out[state, col] = amp
    return out


def independent_fixture_matrix(rungs):
    f = fixture(rungs); n = f['modes']; h = np.zeros((1 << n, 1 << n))
    for term in f['hamiltonian']:
        h += float(Fraction(term['coefficient'])) * car_word_matrix(term['word'], n)
    return h


def original_repulsive_matrix(rungs):
    """Original interleaved-spin repulsive Hubbard Hamiltonian."""
    sites = 2 * rungs; modes = 2 * sites; h = np.zeros((1 << modes, 1 << modes))
    edges = [(2*r, 2*r+1) for r in range(rungs)]
    edges += [(2*r+s, 2*(r+1)+s) for r in range(rungs-1) for s in (0, 1)]
    terms = []
    for a, b in edges:
        for spin in (0, 1):
            for i, j in ((a, b), (b, a)):
                terms.append((-1., ((1, 2*i+spin), (0, 2*j+spin))))
    for i in range(sites):
        terms.append((8., ((1, 2*i), (0, 2*i), (1, 2*i+1), (0, 2*i+1))))
    for c, w in terms: h += c * car_word_matrix(w, modes)
    return h


def particle_hole_signed_permutation(rungs):
    """Map original occupation basis into the transformed basis.

    Down-spin creators become eta_i times down-spin annihilators, acting on
    the reference with every down orbital occupied. The explicit CAR signs are
    retained in the resulting signed permutation.
    """
    sites = 2 * rungs; modes = 2 * sites; n = 1 << modes
    P = np.zeros((n, n)); reference = sum(1 << (2*i+1) for i in range(sites))
    eta = [(-1) ** (i // 2 + i % 2) for i in range(sites)]
    for old in range(n):
        state = reference; phase = 1
        for mode in reversed(range(modes)):
            if not (old & (1 << mode)): continue
            site, spin = divmod(mode, 2)
            create = spin == 0
            op = ((1, mode),) if create else ((0, mode),)
            phase *= eta[site] if spin else 1
            acted = car_word_matrix(op, modes)[:, state]
            if not np.any(acted): phase = 0; break
            state = int(np.flatnonzero(acted)[0]); phase *= int(acted[state])
        if phase: P[state, old] = phase
    return P


def mpo_dense(arrays):
    n = len(arrays); d = arrays[0].shape[2]; out = np.zeros((d ** n, d ** n))
    for bra in range(d ** n):
        bs = [(bra >> i) & 1 for i in range(n)]
        for ket in range(d ** n):
            ks = [(ket >> i) & 1 for i in range(n)]
            env = np.ones((1, 1))
            for w, b, k in zip(arrays, bs, ks): env = env @ w[:, :, b, k]
            out[bra, ket] = env[0, 0]
    return out


class ModelBindingTests(unittest.TestCase):
    def test_transformed_mpo_matches_independent_car_model(self):
        h = independent_fixture_matrix(2)
        arrays, charges = numerical_mpo(mpo(2))
        self.assertEqual(charges[0], [(0, 0)]); self.assertEqual(charges[-1], [(0, 0)])
        self.assertTrue(np.allclose(mpo_dense(arrays), h, atol=1e-12))

    def test_explicit_particle_hole_mapping(self):
        original = original_repulsive_matrix(2)
        transformed = independent_fixture_matrix(2)
        P = particle_hole_signed_permutation(2)
        self.assertTrue(np.allclose(P.T @ P, np.eye(256)))
        # Fixture includes the +8*rungs particle-hole constant.
        sector=[x for x in range(256) if sum((x>>(2*i))&1 for i in range(4))==2
                and sum((x>>(2*i+1))&1 for i in range(4))==2]
        restricted=P @ original @ P.T
        self.assertTrue(np.allclose(restricted[np.ix_(sector,sector)],
                                    transformed[np.ix_(sector,sector)],atol=1e-12))

    def test_trace_seed_norm_and_identity_grouping(self):
        arrays, charges = trace_seed(2)
        # Contract the MPS directly; it represents one identity matrix per
        # conserved (up,down) split, hence the binomial count.
        e = np.array([[1.]])
        for a in arrays: e = np.einsum('ab,asi,bsj->ij', e, a, a)
        self.assertAlmostEqual(float(e[0, 0]), seed_squared_norm(2))
        self.assertEqual(seed_squared_norm(2), 6)
        self.assertEqual(charges[0], [(0, 0)])
        self.assertEqual(charges[-1], [(2, 2)])

    def test_seed_is_signed_identity_after_up_hole_regrouping(self):
        arrays, charges = trace_seed(2)
        # Regroup the interleaved fermions to all-up then all-hole. Each
        # hole preceding an up mode contributes one antisymmetric swap.
        for up in range(16):
            for hole in range(16):
                env=np.ones(1)
                for site in range(4):
                    env=env @ arrays[2*site][:,(up>>site)&1,:]
                    env=env @ arrays[2*site+1][:,(hole>>site)&1,:]
                inversions=sum(((hole>>i)&1)*((up>>j)&1)
                               for i in range(4) for j in range(i+1,4))
                grouped=float(env[0])*(-1)**inversions
                # q=2 gives the single global phase (-1)^(q(q-1)/2)=-1.
                expected=-1. if up==hole and up.bit_count()==2 else 0.
                self.assertAlmostEqual(grouped,expected)



if __name__ == '__main__': unittest.main()
