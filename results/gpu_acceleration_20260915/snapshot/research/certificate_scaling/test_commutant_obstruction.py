"""Independent exact controls for quartic commutant obstruction certificates."""
from __future__ import annotations

import json
import unittest
from fractions import Fraction as F

import numpy as np

from experiments.marginal_symbolic import add, canonical, decode, mono, product, scale, number_shift
from research.certificate_scaling.hidden_density_basis import commutator_column, commutator_map, hidden_fixture, rotate
from research.certificate_scaling.commutant_obstruction import (
    inertia_integer, integer_gram, number_ideal_core, obstruction,
)


def explicit_antisymmetric(poly, i, j):
    q = add(mono(((1, i), (0, j))), scale(mono(((1, j), (0, i))), -1))
    return {w: c for w, c in canonical(add(product(q, poly), scale(product(poly, q), -1))).items() if c}


class CommutantObstructionControls(unittest.TestCase):
    def test_fraction_free_inertia_matches_dense_signs(self):
        cases = [
            [[3, 1, 0], [1, -2, 1], [0, 1, 4]],
            [[2, -1, 1, 0], [-1, 3, 0, 1], [1, 0, -4, 2], [0, 1, 2, 5]],
        ]
        for matrix in cases:
            exact = inertia_integer(matrix)
            values = np.linalg.eigvalsh(np.asarray(matrix, dtype=float))
            self.assertEqual(exact["positive"], int(np.count_nonzero(values > 1e-9)))
            self.assertEqual(exact["negative"], int(np.count_nonzero(values < -1e-9)))

    def test_inertia_on_exact_unit_lower_congruences(self):
        for size in range(2, 8):
            d = [(-1 if (3 * i + size) % 4 == 0 else 1) * (i + 2) for i in range(size)]
            l = [[1 if i == j else (i * 3 - j * 2) if i > j else 0 for j in range(size)]
                 for i in range(size)]
            matrix = [[sum(l[i][k] * d[k] * l[j][k] for k in range(size))
                       for j in range(size)] for i in range(size)]
            got = inertia_integer(matrix)
            self.assertEqual(got["positive"], sum(x > 0 for x in d))
            self.assertEqual(got["negative"], sum(x < 0 for x in d))

    def test_zero_leading_pivot_is_rejected(self):
        with self.assertRaises(ValueError):
            inertia_integer([[0, 1], [1, 0]])

    def test_integer_gram_matches_direct_commutator_gram(self):
        h = hidden_fixture(4)
        gram, denominator, _, _ = integer_gram(h, 4)
        columns, words = commutator_map(h, 4)
        matrix = np.array([[float(col.get(w, 0)) for col in columns] for w in words])
        reconstructed = np.asarray(gram, dtype=float) / denominator**2
        self.assertTrue(np.allclose(reconstructed, matrix.T @ matrix, atol=1e-12))

    def test_number_ideal_core_idempotence_and_density_input(self):
        h = hidden_fixture(4)
        core = number_ideal_core(h, 4)
        self.assertEqual(number_ideal_core(core, 4), core)
        density = product(mono(((1, 0), (0, 0))), mono(((1, 1), (0, 1))))
        density_core = number_ideal_core(density, 4)
        self.assertTrue(all(len(w) == 4 and {w[0][1], w[1][1]} == {w[2][1], w[3][1]}
                            for w in density_core))

    def test_number_ideal_core_removes_lift_and_is_rotation_covariant(self):
        h = hidden_fixture(4)
        q = add(mono(((1, 0), (0, 1))), mono(((1, 1), (0, 0))))
        lift = {w: c for w, c in product(number_shift(4, 0), q).items() if len(w) == 4}
        self.assertEqual(number_ideal_core(add(h, lift), 4), number_ideal_core(h, 4))
        u = [[F(3,5), F(-4,5), F(0), F(0)], [F(4,5), F(3,5), F(0), F(0)],
             [F(0), F(0), F(4,5), F(-3,5)], [F(0), F(0), F(3,5), F(4,5)]]
        ut = [list(row) for row in zip(*u)]
        self.assertEqual(number_ideal_core(rotate(h, u), 4), rotate(number_ideal_core(h, 4), u))

    def test_antisymmetric_commutators_match_direct_car(self):
        h = hidden_fixture(4)
        quartic = {w: c for w, c in h.items() if len(w) == 4}
        for i in range(4):
            for j in range(i + 1, 4):
                self.assertEqual(commutator_column(quartic, i, j, True),
                                 explicit_antisymmetric(quartic, i, j))

    def test_hidden_fixture_has_zero_obstruction_and_molecular_fixture_positive(self):
        hidden = obstruction(hidden_fixture(4), 4)
        self.assertEqual(hidden["frobenius_distance_squared_lower"], "0")
        with open("results/marginal_molecule_stress/h4_square_degree3_certificate.json") as stream:
            source = json.load(stream)
        molecular = obstruction(decode(source["hamiltonian"], source["modes"], 4), source["modes"])
        self.assertEqual(molecular["inertia"]["positive"], 33)
        self.assertEqual(molecular["inertia"]["negative"], 3)
        self.assertEqual(molecular["frobenius_distance_squared_lower"], "1/8000")
        self.assertGreater(molecular["frobenius_distance_lower_float"], 0.0016)

    def test_complex_unitary_obstruction_known_result(self):
        with open("results/marginal_molecule_stress/h4_square_degree3_certificate.json") as stream:
            source = json.load(stream)
        h = decode(source["hamiltonian"], source["modes"], 4)
        result = obstruction(h, source["modes"], complex_unitary=True)
        self.assertEqual(result["total_negative_inertia"], 4)
        self.assertEqual(result["frobenius_distance_squared_lower"], "1/10000")
        self.assertAlmostEqual(result["frobenius_distance_lower_float"], 0.01, places=12)

    def test_complex_quotient_bound_and_known_approximation(self):
        from experiments.marginal_symbolic import adj
        w=mono(((1,0),(1,1),(0,2),(0,3)),F(1,10))
        perturbation=canonical(add(w,adj(w)))
        h=add(hidden_fixture(4),perturbation)
        known_squared=sum(c*c for c in perturbation.values())
        result=obstruction(h,4,F(1,1000000),True,True)
        self.assertLessEqual(F(result['frobenius_distance_squared_lower']),known_squared)
        with open('results/certificate_scaling/active_space_ladder/h4/fixture.json') as stream:
            source=json.load(stream)
        molecular=obstruction(decode(source['hamiltonian'],8,4),8,complex_unitary=True,quotient_number_ideal=True)
        self.assertEqual(molecular['total_negative_inertia'],4)
        self.assertEqual(molecular['frobenius_distance_squared_lower'],'1/10000')

    def test_real_input_gate_rejects_nonhermitian_projection(self):
        bad=add(hidden_fixture(4),mono(((1,0),(0,1)),F(1,7)))
        with self.assertRaises(ValueError):number_ideal_core(bad,4)
        with self.assertRaises(ValueError):obstruction(bad,4,complex_unitary=True)


if __name__ == "__main__":
    unittest.main()
