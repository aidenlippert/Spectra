import unittest
from fractions import Fraction as F
import numpy as np
from .local_basis_state import local_fock_gate,apply_gate,givens_design
from .orbital_order import swap_neighbors
from .validated_tensor import right_canonicalize,tensor_bounds
from .test_validated_tensor import dense


class TestCircuit(unittest.TestCase):
 def test_exact_local_unitarity_and_generator(self):
  g=local_fock_gate(F(3,5),F(4,5))
  for i in range(16):
   for j in range(16):self.assertEqual(sum(g[k][i]*g[k][j] for k in range(16)),int(i==j))
  # one-alpha-electron block: |1000>, |0010>
  self.assertEqual(g[8][8],F(3,5));self.assertEqual(g[2][8],F(-4,5))
  self.assertEqual(g[15][15],F(1))
 def test_direct_gate_and_permutation(self):
  rng=np.random.default_rng(190)
  ts=[rng.normal(size=(1,4,3))+1j*rng.normal(size=(1,4,3)),rng.normal(size=(3,4,1))+1j*rng.normal(size=(3,4,1))]
  g=local_fock_gate(F(3,5),F(4,5));out,e=apply_gate(ts,0,g)
  expected=np.array(g,float)@dense(ts)
  self.assertLessEqual(np.linalg.norm(dense(out)-expected),e)
  out,e=swap_neighbors(ts,0);old=dense(ts).reshape(4,4);expected=old.T.copy()
  for i in (1,2):
   for j in (1,2):expected[i,j]*=-1
  self.assertLessEqual(np.linalg.norm(dense(out)-expected.reshape(-1)),e)
 def test_checked_canonicalization(self):
  rng=np.random.default_rng(21)
  ts=[rng.normal(size=(1,2,3)),rng.normal(size=(3,2,1))]
  out,e=right_canonicalize(ts)
  self.assertLessEqual(np.linalg.norm(dense(ts)-dense(out)),e)
  replay,er=right_canonicalize(ts,proposal=out)
  self.assertLessEqual(np.linalg.norm(dense(ts)-dense(replay)),er)

if __name__=='__main__':unittest.main()
