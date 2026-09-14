import copy
from fractions import Fraction as F
from itertools import combinations
import unittest
from experiments.marginal_hubbard_polynomial import build
from experiments.marginal_symmetry_moments import SymmetryMomentOracle
from experiments.marginal_determinant_tree import DeterminantOracle
from experiments.marginal_symbolic import decode,encode,add,mono


def states():
    return [sum(1<<(2*i) for i in a)+sum(1<<(2*i+1) for i in b)
            for a in combinations(range(4),2) for b in combinations(range(4),2)]
def act(o,v):
    out={}
    for s,a in v.items():
        for t,b in o.action(s).items():out[t]=out.get(t,0)+a*b
    return {s:a for s,a in out.items() if a}
def expand(o,v):
    out={}
    for r,a in v.items():
        for s,p in o.orbit(r)[3].items():out[s]=out.get(s,0)+a*p
    return {s:a for s,a in out.items() if a}

class SymmetryMomentsTests(unittest.TestCase):
    def test_all_four_site_signed_orbit_actions(self):
        h=build(4,4,1);o=SymmetryMomentOracle(h);bare=DeterminantOracle(h);reps=set()
        for s in states():
            f,pf=o.flip(s);ff,pff=o.flip(f);self.assertEqual((ff,pf*pff),(s,1))
            c,pc=o.particle_hole(s);cc,pcc=o.particle_hole(c);self.assertEqual((cc,pc*pcc),(s,1))
            fc,pfc=o.flip(c);cf,pcf=o.particle_hole(f);self.assertEqual((fc,pc*pfc),(cf,pf*pcf))
            entry=o.orbit(s)
            if entry is not None:reps.add(entry[0])
        self.assertLess(len(reps),len(states()))
        for r in reps:
            v=expand(o,{r:1})
            self.assertEqual(expand(o,o.action(r)),act(bare,v))
            self.assertEqual(o.compress(v),{r:1})
        for r in reps:
            for s in reps:self.assertEqual(o.orbit(r)[2]*o.action(s).get(r,0),o.orbit(s)[2]*o.action(r).get(s,0))
    def test_moments_through_fourteen_against_direct_car(self):
        h=build(4,4,1);o=SymmetryMomentOracle(h);bare=DeterminantOracle(h)
        reps=sorted({o.orbit(s)[0] for s in states() if o.orbit(s) is not None})[:3]
        V=[expand(o,{r:1}) for r in reps];moments,receipt=o.moments(V,14)
        power=V
        for k in range(15):
            actual=[[sum(a*v.get(s,0) for s,a in u.items()) for v in power] for u in V]
            self.assertEqual(moments[k],actual)
            if k<14:power=[act(bare,v) for v in power]
        self.assertLess(receipt['unique_determinant_source_actions'],len(bare.cache))
    def test_bad_symmetry_and_boundary_refused(self):
        h=build(4,4,1);o=SymmetryMomentOracle(h)
        r=next(s for s in states() if o.orbit(s) is not None and o.orbit(s)[2]>1)
        with self.assertRaisesRegex(ValueError,'not invariant'):o.compress({r:1})
        with self.assertRaisesRegex(ValueError,'exact integer'):o.compress({r:1.0})
        with self.assertRaises(ValueError):o.moments([],2)
        with self.assertRaises(ValueError):o.moments([expand(o,{o.orbit(r)[0]:1})],25)
        bad=copy.deepcopy(h);bad['particles']=2
        with self.assertRaisesRegex(ValueError,'half filling'):SymmetryMomentOracle(bad)
        bad=copy.deepcopy(h);poly=decode(bad['hamiltonian'],8,4)
        bad['hamiltonian']=encode(add(poly,mono(((1,0),(0,0)))))
        with self.assertRaisesRegex(ValueError,'spin-exchange'):SymmetryMomentOracle(bad)
        bad['hamiltonian']=encode(add(poly,mono(((1,0),(0,0))),mono(((1,1),(0,1)))))
        with self.assertRaisesRegex(ValueError,'Particle-hole'):SymmetryMomentOracle(bad)
    def test_nonuniform_bipartite_ring_moments(self):
        h=build(4,4,1);poly=decode(h['hamiltonian'],8,4)
        extra=[mono(((1,s),(0,6+s)),F(-2,3)) for s in (0,1)]
        extra += [mono(((1,6+s),(0,s)),F(-2,3)) for s in (0,1)]
        h['hamiltonian']=encode(add(poly,*extra));o=SymmetryMomentOracle(h);bare=DeterminantOracle(h)
        reps=sorted({o.orbit(s)[0] for s in states() if o.orbit(s) is not None})[:2]
        V=[expand(o,{r:1}) for r in reps];moments,_=o.moments(V,8);power=V
        for k in range(9):
            actual=[[sum(a*v.get(s,0) for s,a in u.items()) for v in power] for u in V]
            self.assertEqual(moments[k],actual)
            if k<8:power=[act(bare,v) for v in power]
if __name__=='__main__':unittest.main()
