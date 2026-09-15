from fractions import Fraction
import unittest
from experiments.pauli import multiply, commutator_i, closure, generator_matrix, commutation_adjacency

class PauliTests(unittest.TestCase):
  def test_exact_products_and_commutator(self):
    assert multiply("X", "Y") == (1j, "Z")
    self.assertEqual(commutator_i({"X": Fraction(1, 3)}, "Y"), {"Z": Fraction(-2, 3)})
    self.assertEqual(commutator_i({"X": 2**80}, "Y"), {"Z": Fraction(-(2**81))})

  def test_closure_and_budget(self):
    b, complete = closure([{"X": 1}], ["Y"], budget=8)
    self.assertTrue(complete); self.assertEqual(b, {"Y", "Z"})
    b, complete = closure([{"X": 1}], ["Y"], budget=1)
    self.assertFalse(complete); self.assertEqual(b, {"Y"})

  def test_residual_and_adjacency(self):
    labels, matrix, residual = generator_matrix({"X": 1, "Z": 1}, ["Y"])
    self.assertEqual(labels, ["Y"]); self.assertEqual(matrix, [[0]]); self.assertEqual(residual["Y"], 4)
    self.assertEqual(commutation_adjacency(["X", "Y", "Z"]), 3)

  def test_rejects_bad_inputs(self):
    for bad in [("X", "YY"), ("A", "X")]:
      with self.assertRaises(ValueError): multiply(*bad)

  def test_row_orientation(self):
    labels, matrix, _ = generator_matrix({"X": 1}, ["Y", "Z"])
    self.assertEqual(labels, ["Y", "Z"])
    self.assertEqual(matrix, [[0, -2], [2, 0]])

if __name__ == "__main__": unittest.main()
