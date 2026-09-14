"""Dependency-free tests of symbolic certificate soundness and rejection paths."""
from fractions import Fraction as F
import json
from pathlib import Path
import unittest

from experiments.marginal_symbolic import (
    encode, expand_orbits, model_hamiltonian, number_shift, verify, verified_residual, expand_squares,
)
from experiments.marginal_hunt_car import add, mono


def number_certificate(modes, particles):
    return {"modes": modes, "particles": particles,
            "hamiltonian": encode(number_shift(modes, 0)),
            "number_multiplier": encode(mono(())), "b": str(particles),
            "denominator": 1, "blocks": []}


def triple_certificate():
    blocks = []
    for triple in ((0, 1, 2), (3, 4, 5)):
        word = tuple((0, i) for i in triple)
        for w in (word, tuple((1-c, i) for c, i in reversed(word))):
            blocks.append({"words": [w], "factor": [[1]]})
    return {"modes": 6, "particles": 3, "hamiltonian": encode(model_hamiltonian(0)),
            "number_multiplier": encode(mono(())), "b": "1", "denominator": 1,
            "blocks": blocks}


class SymbolicTests(unittest.TestCase):
    def test_gram_symmetry_expansion_matches_individual_squares(self):
        from experiments.marginal_symbolic import product, scale, canonical
        from experiments.marginal_hunt_car import adj
        words=[((0,0),),((0,1),),((1,2),(0,1),(0,0))]
        factors=[[2,-3,1],[0,4,-2],[-5,1,3]]
        blocks=[{'words':words,'factor':factors}]
        expected={}
        for row in factors:
            q=add(*(scale(mono(w),F(c,7)) for w,c in zip(words,row)))
            expected=add(expected,product(canonical(adj(q)),q))
        actual,stats=expand_squares(blocks,7,3)
        self.assertEqual(actual,expected)
        self.assertEqual(stats['factor_rows'],3)
        self.assertEqual(expand_squares([{'words':words,'factor':[]}],7,3)[0],{})

    def test_shared_verified_residual_preserves_receipt(self):
        cert=triple_certificate();residual,receipt=verified_residual(cert)
        self.assertEqual(residual,{})
        self.assertEqual(receipt,verify(cert))

    def test_number_identity_including_large_and_edge_sectors(self):
        for modes, particles in ((1, 0), (1, 1), (6, 3), (30, 15), (1000, 500)):
            result = verify(number_certificate(modes, particles))
            self.assertEqual(result["residual_l1"], "0")
            self.assertEqual(result["lower"], str(particles))

    def test_exact_triple_identity(self):
        result = verify(triple_certificate())
        self.assertEqual(result["lower"], "1")
        self.assertEqual(result["residual_l1"], "0")

    def test_changed_number_is_not_silently_accepted(self):
        c = triple_certificate()
        c["particles"] = 2
        result = verify(c)
        self.assertEqual(result["lower"], "0")
        self.assertEqual(result["residual_l1"], "1")

    def test_scalar_and_factor_corruption_charged_to_residual(self):
        c = triple_certificate()
        c["b"] = "1001"
        self.assertEqual(verify(c)["lower"], "1")
        c = triple_certificate()
        c["blocks"][0]["factor"] = [[2]]
        result = verify(c)
        self.assertGreater(F(result["residual_l1"]), 0)
        self.assertLessEqual(F(result["lower"]), 1)

    def test_invalid_coefficients_and_nonhermitian_multiplier_rejected(self):
        c = triple_certificate()
        c["number_multiplier"] = encode(mono(((1, 0), (0, 1))))
        with self.assertRaises(ValueError):
            verify(c)
        c["number_multiplier"] = encode(add(mono(((1, 0),)), mono(((0, 0),))))
        with self.assertRaises(ValueError):
            verify(c)
        c = triple_certificate()
        c["blocks"][0]["factor"] = [[.5]]
        with self.assertRaises(ValueError):
            verify(c)
        c = triple_certificate()
        c["number_multiplier"][0]["coefficient"] = .1
        with self.assertRaises(ValueError):
            verify(c)

    def test_saved_hopping_certificate_without_numpy(self):
        root = Path(__file__).resolve().parents[1]
        c = json.loads((root / "results/marginal_symbolic/full_mixed_matched.json").read_text())
        result = verify(c)
        self.assertGreater(F(result["lower"]), F(5509, 10000))
        self.assertLess(F(result["residual_l1"]), F(1, 10000))

    def test_compressed_multiplier_and_malformed_orbit(self):
        root = Path(__file__).resolve().parents[1]
        c = json.loads((root / "results/marginal_symbolic/matched_symmetry_compressed.json").read_text())
        self.assertEqual(len(c["multiplier_orbits"]["orbits"]), 10)
        self.assertEqual(len(expand_orbits(c["multiplier_orbits"], 6)), 64)
        result = verify(c)
        self.assertGreater(F(result["lower"]), F(5509, 10000))
        c["multiplier_orbits"]["permutations"][0] = [0]*6
        with self.assertRaises(ValueError):
            verify(c)


if __name__ == "__main__":
    unittest.main()
