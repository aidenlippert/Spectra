from copy import deepcopy
from fractions import Fraction as F
from itertools import product
import unittest
import numpy as np
from experiments.v3_locality import (certify, verify, chain_hamiltonian, neighborhood,
                                     exp_negative_interval, propagate_enclosed, reduced_columns)
from experiments.pauli import commutator_i


class LocalityTests(unittest.TestCase):
    def test_exact_physical_certificate_and_rejections(self):
        h, seed = chain_hamiltonian(6)
        c = certify(h, seed)
        self.assertTrue(verify(h, c))
        self.assertLess(F(c['uniform_all_time_error_bound']), F(1, 1000))
        for field, value in [('depth', 3.5), ('kappa', '0'), ('uniform_all_time_error_bound', '0')]:
            bad = deepcopy(c); bad[field] = value
            self.assertFalse(verify(h, bad))
        bad = deepcopy(c); bad['basis'].remove(seed)
        self.assertFalse(verify(h, bad))
        self.assertEqual(certify(h, seed, gamma=F(1))['status'], 'no_certificate')
        self.assertEqual(certify(h, seed, cap=1)['status'], 'no_certificate')

    def test_local_certificate_size_stabilizes(self):
        sizes = []
        for n in (12, 24, 48):
            h, seed = chain_hamiltonian(n)
            c = certify(h, seed)
            self.assertTrue(verify(h, c))
            sizes.append(len(c['basis']))
        self.assertEqual(len(set(sizes)), 1)

    def test_exponential_interval_and_reduced_enclosure(self):
        for x in (F(0), F(1, 2), F(10)):
            low, high, _ = exp_negative_interval(x)
            self.assertLessEqual(float(low) - 1e-15, np.exp(-float(x)))
            self.assertGreaterEqual(float(high) + 1e-15, np.exp(-float(x)))
        h, seed = chain_hamiltonian(3)
        c = certify(h, seed)
        result = propagate_enclosed(h, c, F(1, 1000))
        basis = result['basis']; columns = reduced_columns(h, basis, F(c['gamma']))
        a = np.zeros((len(basis), len(basis)))
        for j, col in enumerate(columns):
            for i, value in col.items(): a[i, j] = float(value)
        eig, vectors = np.linalg.eig(a.astype(complex))
        initial = np.array([int(w == seed) for w in basis])
        dense = vectors @ (np.exp(eig / 1000) * np.linalg.solve(vectors, initial))
        center = np.array([float(v) for v in result['center']])
        self.assertLess(np.linalg.norm(dense - center, 1), float(result['numerical_error_bound']) + 1e-12)

    def test_interactions_cannot_be_ignored_at_requested_accuracy(self):
        h, seed = chain_hamiltonian(3)
        c = certify(h, seed)
        pred = propagate_enclosed(h, c, F(1, 1000))
        axes = {0: 'X', 1: 'Y', 2: 'X'}
        value = sum(v for w, v in zip(pred['basis'], pred['center'])
                    if all(p == 'I' or axes[i] == p for i, p in enumerate(w)))
        lower = value - F(c['uniform_all_time_error_bound']) - pred['numerical_error_bound']
        # The same physical state has zero expectation under damping-only motion.
        self.assertGreater(lower, F(c['tolerance']))

    def test_full_pauli_reference_obeys_bound(self):
        h, seed = chain_hamiltonian(3)
        c = certify(h, seed)
        full = [''.join(x) for x in product('IXYZ', repeat=3) if ''.join(x) != 'III']
        columns = reduced_columns(h, full, F(c['gamma']))
        a = np.zeros((len(full), len(full)))
        for j, col in enumerate(columns):
            for i, value in col.items(): a[i, j] = float(value)
        eig, vec = np.linalg.eig(a.astype(complex))
        initial = np.array([int(w == seed) for w in full])
        eig_initial = np.linalg.solve(vec, initial)
        for t in (F(1, 1000), F(1, 100), F(1, 32)):
            ref = vec @ (np.exp(eig * float(t)) * eig_initial)
            approx = propagate_enclosed(h, c, t)
            embedded = np.zeros(len(full))
            for word, value in zip(approx['basis'], approx['center']): embedded[full.index(word)] = float(value)
            error = np.linalg.norm(ref - embedded, 1)
            self.assertLessEqual(error, float(F(c['uniform_all_time_error_bound']) + approx['numerical_error_bound']) + 1e-12)


if __name__ == '__main__': unittest.main()
