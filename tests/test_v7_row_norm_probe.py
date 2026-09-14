import unittest
from fractions import Fraction as F
from experiments.v7_row_norm_probe import row_bound
class RowNormTests(unittest.TestCase):
 def test_phase_cancellation_and_complex_entries(self):
  self.assertEqual(row_bound({'XX':F(1),'YY':F(1)},2)[0],2)
  self.assertEqual(row_bound({'XY':F(1),'YX':F(1)},2)[0],2)
  b,c=row_bound({'X':F(1),'Y':F(1)},1)
  self.assertGreaterEqual(b*b,F(2))
  self.assertLess(b,F(14142135623730952,10**16))
  self.assertEqual(c['pauli_entry_updates'],4)
 def test_exponential_budget_refuses(self):
  with self.assertRaises(ValueError):row_bound({'X'*7:F(1)},7)
if __name__=='__main__':unittest.main()
