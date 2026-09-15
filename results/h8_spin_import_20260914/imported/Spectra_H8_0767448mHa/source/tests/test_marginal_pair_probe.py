from fractions import Fraction as F
from itertools import combinations
import unittest

from experiments.marginal_pair_probe import pair_polynomials, probe_block, solve
from experiments.marginal_symbolic import add, adj, encode, mono, number_shift, product, scale, verify
from tests.test_marginal_hunt_car import act


class PairProbeTests(unittest.TestCase):
    def test_span_discovery_exports_valid_quartic_certificate(self):
        cert,receipt=solve(modes=4,kind='span',solver='CLARABEL')
        self.assertEqual(cert['operator_degree'],4)
        self.assertLess(F(receipt['width']),F(1,1000))

    def test_requested_square_is_diagonal_not_pair_coherence(self):
        for a in pair_polynomials(6):
            self.assertEqual(product(a,a),{})
            b=add(a,scale(adj(a),-1))
            square=product(adj(b),b)
            self.assertEqual(square,add(product(adj(a),a),product(a,adj(a))))
            for state in range(1<<6):
                self.assertTrue(set(act(square,state)).issubset({state}))

    def test_polynomial_transform_reproduces_requested_operators(self):
        block,transform=probe_block(6,'antisymmetric')
        expected=[add(a,scale(adj(a),-1)) for a in pair_polynomials(6)]
        for j,p in enumerate(expected):
            actual={w:F(int(transform[i,j])) for i,w in enumerate(block['words']) if transform[i,j]}
            self.assertEqual(actual,p)

    def test_quartic_square_requires_opt_in_and_retains_degree_eight(self):
        block,tr=probe_block(4,'antisymmetric')
        block['factor']=[[int(x) for x in tr[:,0]]]
        cert={'modes':4,'particles':2,'hamiltonian':[], 'number_multiplier':[],
              'b':'0','denominator':1,'blocks':[block]}
        with self.assertRaises(ValueError):verify(cert)
        cert['operator_degree']=4
        receipt=verify(cert)
        self.assertEqual(receipt['residual_max_degree'],8)
        self.assertLess(F(receipt['lower']),0)
        cert['operator_degree']=5
        with self.assertRaises(ValueError):verify(cert)

    def test_degree_six_number_multiplier_cancels_quartic_factor(self):
        shift=number_shift(6,3)
        # Y=n0+n1 has Y^2 of degree four, so X=-(Nhat-N)Y^2
        # genuinely has degree six. B=(Nhat-N)Y has degree four.
        y=add(mono(((1,0),(0,0))),mono(((1,1),(0,1))))
        b=product(shift,y)
        x=scale(product(shift,product(y,y)),-1)
        self.assertEqual(max(map(len,x)),6)
        words=list(b)
        cert={'modes':6,'particles':3,'hamiltonian':[], 'number_multiplier':encode(x),
              'b':'0','denominator':1,'operator_degree':4,
              'blocks':[{'words':words,'factor':[[int(b[w]) for w in words]]}]}
        self.assertEqual(F(verify(cert)['lower']),0)
        for indices in combinations(range(6),3):
            self.assertEqual(act(b,sum(1<<i for i in indices)),{})


if __name__=='__main__':unittest.main()
