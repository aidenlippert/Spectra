"""Independent dense controls for the restricted stoquastic chain certificate."""
from __future__ import annotations

import itertools
import unittest
from fractions import Fraction as F

import numpy as np

from research.certificate_scaling.stoquastic_chain_certificate import (
    index, local_energy, minmax_dp, parent_hamiltonian, verify,
)


def dense_hamiltonian(h):
    n = h["sites"]
    dim = 1 << n
    mat = np.zeros((dim, dim))
    for state in range(dim):
        bits = [(state >> i) & 1 for i in range(n)]
        for i, table in enumerate(h["diagonal_tables"]):
            mat[state, state] += float(F(table[index(bits[(i - 1) % n], bits[i], bits[(i + 1) % n])]))
            flipped = state ^ (1 << i)
            mat[flipped, state] += -float(F(h["flip_rates"][i]))
    return mat


class StoquasticChainControls(unittest.TestCase):
    def test_dp_matches_bruteforce_and_dense_spectrum_is_bracketed(self):
        for n in (4, 5):
            h = parent_hamiltonian(n)
            for q in (F(2), F(3, 2), F(5, 3)):
                tables = local_energy(h, q)
                lo, hi, _ = minmax_dp(tables)
                brute = [sum(tables[i][index(x[(i - 1) % n], x[i], x[(i + 1) % n])]
                              for i in range(n))
                         for x in itertools.product(range(2), repeat=n)]
                self.assertEqual((lo, hi), (min(brute), max(brute)))
                rec = verify({"hamiltonian": h,
                              "amplitude_domain_wall_weight": str(q)})
                eig = np.linalg.eigvalsh(dense_hamiltonian(h))
                self.assertLessEqual(float(F(rec["lower"])), eig[0] + 1e-9)
                self.assertGreaterEqual(float(F(rec["upper"])), eig[0] - 1e-9)

    def test_perturbed_diagonal_and_rates_remain_bracketed(self):
        h = parent_hamiltonian(5)
        h["diagonal_tables"][2][index(0, 1, 0)] = "7/5"
        h["flip_rates"] = ["1", "3/2", "1/2", "5/4", "2"]
        eig = np.linalg.eigvalsh(dense_hamiltonian(h))
        rec = verify({"hamiltonian": h, "amplitude_domain_wall_weight": "3/2"})
        self.assertLessEqual(float(F(rec["lower"])), eig[0] + 1e-9)
        self.assertGreaterEqual(float(F(rec["upper"])), eig[0] - 1e-9)

    def test_local_terms_are_noncommuting(self):
        h = parent_hamiltonian(4)
        # Independent local-term construction: diagonal term at site i and
        # flip term at i. Adjacent terms share a site and do not commute.
        n = 4; dim = 1 << n
        local = []
        for i, table in enumerate(h["diagonal_tables"]):
            d = np.zeros((dim, dim)); f = np.zeros((dim, dim))
            for state in range(dim):
                b = [(state >> j) & 1 for j in range(n)]
                d[state, state] = float(F(table[index(b[(i-1)%n], b[i], b[(i+1)%n])]))
                f[state ^ (1 << i), state] = -1
            local.append((d, f))
        self.assertGreater(np.linalg.norm(local[0][0] @ local[1][1] - local[1][1] @ local[0][0]), 0)

    def test_nonstoquastic_rate_rejected(self):
        h = parent_hamiltonian(4); h["flip_rates"][0] = "-1"
        with self.assertRaises(ValueError):
            verify({"hamiltonian": h, "amplitude_domain_wall_weight": "2"})


if __name__ == "__main__":
    unittest.main()
