"""Independent CAR and transport controls for hidden density-basis discovery."""
from __future__ import annotations

import unittest
from fractions import Fraction as F

from experiments.marginal_symbolic import add, canonical, mono, product, scale, encode
from research.certificate_scaling.hidden_density_basis import (
    commutator_column, generators, hidden_fixture, rotate, discover, replay, detect,
)
from research.certificate_scaling.fermionic_ratio_chain import hamiltonian


def explicit_commutator(poly, i, j):
    q = {}
    for p in range(4):
        if p == i == j:
            q = add(q, mono(((1, p), (0, p))))
        elif p in (i, j):
            q = add(q, mono(((1, p), (0, (j if p == i else i)))))
    return {w: c for w, c in canonical(add(product(q, poly), scale(product(poly, q), -1))).items() if c}


def explicit_rotate(poly, u):
    out = {}
    for word, coefficient in poly.items():
        terms = {(): coefficient}
        for cr, old in word:
            expanded = {}
            for prefix, value in terms.items():
                for p in range(len(u)):
                    letter = (cr, p)
                    expanded[ prefix + (letter,) ] = expanded.get(prefix + (letter,), F(0)) + value * u[p][old]
            terms = expanded
        for word2, value in terms.items():
            for canonical_word, sign in canonical(mono(word2)).items():
                out[canonical_word] = out.get(canonical_word, F(0)) + value * sign
    return {w: c for w, c in out.items() if c}


class HiddenDensityControls(unittest.TestCase):
    def test_commutator_derivation_all_generators(self):
        h = hidden_fixture(4)
        quartic = {w: c for w, c in h.items() if len(w) == 4}
        for i, j in generators(4):
            self.assertEqual(commutator_column(quartic, i, j),
                             explicit_commutator(quartic, i, j))

    def test_rotate_matches_direct_substitution_and_inverse(self):
        h = hamiltonian(4, F(4))
        # Two rational Givens rotations catch mixed quadratic terms and
        # exterior-square cross terms that a permutation would miss.
        u = [[F(3,5), F(-4,5), F(0), F(0)], [F(4,5), F(3,5), F(0), F(0)],
             [F(0), F(0), F(4,5), F(-3,5)], [F(0), F(0), F(3,5), F(4,5)]]
        self.assertEqual(rotate(h, u), explicit_rotate(h, u))
        transpose_u = [list(row) for row in zip(*u)]
        self.assertEqual(rotate(rotate(h, u), transpose_u), h)

    def test_m4_commutant_modular_rank_gives_nullity_upper_bound(self):
        _, receipt = detect(hidden_fixture(4), 4)
        self.assertEqual(receipt["exact_rank_lower_bound"], 6)
        self.assertEqual(receipt["exact_nullity_upper_bound"], 4)
        self.assertEqual(receipt["numeric_nullity"], 4)

    def test_hidden_m4_discover_and_exact_transport_replay(self):
        h = hidden_fixture(4)
        cert, rec = discover(h, 4, 2)
        self.assertTrue(rec["exact_orthogonality"])
        self.assertTrue(rec["exact_hamiltonian_transport"])
        self.assertEqual(replay(cert)["width"], "0")

    def test_rejects_bad_rotation_sector_hamiltonian_and_inner_upper(self):
        h = hidden_fixture(4); cert, _ = discover(h, 4, 2)
        bad = dict(cert); bad["rotation"] = [row[:] for row in cert["rotation"]]
        bad["rotation"][0][0] = str(F(bad["rotation"][0][0]) + 1)
        with self.assertRaises(ValueError): replay(bad)
        bad = dict(cert); bad["modes"] = 6
        with self.assertRaises(ValueError): replay(bad)
        bad = dict(cert); bad["hamiltonian"] = encode(hamiltonian(4, F(4)))
        with self.assertRaises(ValueError): replay(bad)
        bad = dict(cert); bad["inner_certificate"] = dict(cert["inner_certificate"])
        bad["inner_certificate"]["structured_upper"] = {"kind": "positive_fixed_number_site_ratio_v1", "r": "4"}
        with self.assertRaises(ValueError): replay(bad)


if __name__ == "__main__":
    unittest.main()
