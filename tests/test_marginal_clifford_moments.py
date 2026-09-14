import unittest
from fractions import Fraction as F
from experiments.marginal_hubbard_polynomial import build
from experiments.marginal_clifford_moments import CliffordMomentOracle
from experiments.marginal_symmetry_moments import SymmetryMomentOracle
from experiments.marginal_symbolic import add,mono,adj,encode,decode,scale
from experiments.marginal_transfer_verify import apply_word


def physical_action(h,v):
    out={}
    for s,a in v.items():
        for word,b in h.items():
            z=apply_word(word,s)
            if z:out[z[0]]=out.get(z[0],0)+a*b*z[1]
    return {s:a for s,a in out.items() if a}
def dimer_boundary():return {0x99:1,0x96:-1,0x69:-1,0x66:1}
def general_model():
    h=build(4,4,1);base=decode(h['hamiltonian'],8,4)
    q=mono(((1,0),(1,3),(0,5),(0,6)),F(1,7));flip=mono(((1,0),(0,1)),F(1,5))
    h['hamiltonian']=encode(add(base,q,adj(q),flip,adj(flip),mono(((1,0),(0,0)),F(2,3))))
    return h

class CliffordMomentsTests(unittest.TestCase):
    def test_complete_four_site_frame_action_equals_car(self):
        # Every frame coordinate, including those that mix number sectors.
        for h in [build(4,4,1),general_model()]:
            o=CliffordMomentOracle(h)
            for s in range(256):
                physical=o._change_basis({s:1},False)
                self.assertEqual(o._change_basis(physical,True),{s:1})
                actual=o._change_basis(physical_action(o.base.h,physical),True)
                self.assertEqual(o.action(s),actual)
    def test_dimer_moments_and_general_hamiltonian_transfer(self):
        for h in [build(4,4,1),general_model()]:
            o=CliffordMomentOracle(h);v=dimer_boundary()
            self.assertEqual(o.to_frame(v),{0:2})
            M,r=o.dimer_moments(14);explicit,_=o.moments([v],14);self.assertEqual(M,explicit)
            power=v
            for k in range(15):
                self.assertEqual(M[k][0][0],sum(a*power.get(s,0) for s,a in v.items()))
                if k<14:power=physical_action(o.base.h,power)
            self.assertEqual(r['initial_frame_supports'],[1])
        with self.assertRaises(ValueError):SymmetryMomentOracle(general_model())
    def test_refusal_gates(self):
        o=CliffordMomentOracle(build(4,4,1))
        with self.assertRaises(ValueError):o.to_frame({0:1})
        with self.assertRaisesRegex(ValueError,'exact'):o.to_frame({0x99:1.0})
        with self.assertRaises(ValueError):o.dimer_moments(25)
        with self.assertRaises(ValueError):o.action(256)
        o=CliffordMomentOracle(build(8,4,1))
        for s in range(4096):o.action(s)
        with self.assertRaisesRegex(ValueError,'action budget'):o.action(4096)
if __name__=='__main__':unittest.main()
