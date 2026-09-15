import copy
from fractions import Fraction as F
import unittest

from experiments.marginal_symmetry_transfer import hopping_perturbation,integer_upper,replay
from experiments.marginal_symbolic import decode,encode,add,mono,scale,product,adj


class SymmetryTransferTests(unittest.TestCase):
    def test_hopping_transfer_identity_is_exact(self):
        for modes in (4,6,10):
            for epsilon in (F(0),F(1,1000),F(-1,7)):
                delta,squares,shift=hopping_perturbation(modes,epsilon)
                rhs=mono((),-shift) if shift else {}
                for item in squares:
                    p=decode(item['polynomial'],modes,4)
                    rhs=add(rhs,scale(product(adj(p),p),F(item['weight'])))
                self.assertEqual(rhs,delta)

    def certificate(self):
        h={((1,0),(0,1)):F(-1),((1,1),(0,0)):F(-1)}
        return {'modes':2,'particles':1,'operator_degree':3,'hamiltonian':encode(h),
                'number_multiplier':[],'b':'0','permutations':[[0,1]],'orbit_squares':[],
                'independent_upper':integer_upper(h,2,1)}

    def test_upper_is_exactly_recomputed_not_trusted(self):
        c=self.certificate();c['independent_upper']['numerical_energy']=-999
        c['independent_upper']['norm']=0;c['independent_upper']['upper']=-999
        self.assertEqual(F(replay(c)['upper']),-1)
        c['hamiltonian']=encode({((1,0),(0,1)):F(1),((1,1),(0,0)):F(1)})
        self.assertEqual(F(replay(c)['upper']),1)

    def test_upper_rejects_invalid_vectors(self):
        c=self.certificate()
        for amps in ([0,0],[1],[True,1],[1.5,2]):
            bad=copy.deepcopy(c);bad['independent_upper']['amplitudes']=amps
            with self.assertRaises(ValueError):replay(bad)


if __name__=='__main__':unittest.main()
