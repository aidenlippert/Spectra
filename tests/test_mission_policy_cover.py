import unittest
from fractions import Fraction as F
from experiments.mission_policy_cover import Sample,certify,construct,box_for


class PolicyCoverTests(unittest.TestCase):
    def test_holes_cannot_prove_impossibility(self):
        s=[Sample((F(0),),((F(0),F(0)),(F(1),F(1))))]
        with self.assertRaisesRegex(ValueError,'incomplete'):certify(['0'],s,[1,0],1,F(1,10))
        with self.assertRaisesRegex(ValueError,'overlapping'):certify(['','0'],s,[1,0],1,F(1,10))
        self.assertEqual(certify([''],s,[1,0],1,F(1,10))['status'],'no_feasible_policy_in_covered_domain')

    def test_unobserved_region_remains_possible(self):
        # g(x)=1-2x is consistent: observing failure at x=0 cannot rule out x=1.
        s=[Sample((F(0),),((F(0),F(0)),(F(1),F(1))))]
        self.assertEqual(certify([''],s,[0,2],1,F(1,10))['status'],'unresolved')
        s.append(Sample((F(1),),((F(0),F(0)),(F(-1),F(-1)))))
        r=certify([''],s,[0,2],1,0)
        self.assertEqual(r['status'],'conditional_certificate')
        self.assertEqual(r['policy'],(F(1),))

    def test_point_optimum_is_not_global_lower_bound(self):
        # Known J(x)=x; only x=1 observed. A singleton observation isn't a cover.
        s=[Sample((F(1),),((F(1),F(1)),))]
        r=certify([''],s,[1],1,F(1,10))
        self.assertEqual(r['status'],'unresolved');self.assertEqual(r['lower'],0)
        self.assertEqual(r['upper'],1)

    def test_contradictory_evidence_and_noise_floor(self):
        s=[Sample((F(0),),((F(0),F(0)),)),Sample((F(1),),((F(2),F(2)),))]
        self.assertEqual(certify([''],s,[1],1,F(1))['status'],'inconsistent_assumptions')
        r=construct(lambda p,w:((F(0),F(0)),(F(-1,10),F(1,10))),[0,0],1,F(1,10),8)
        self.assertEqual(r['receipt']['status'],'unresolved');self.assertEqual(r['observation_calls'],8)

    def test_exact_construction_and_replay(self):
        # Mathematical verification of continuum search, not a materials result.
        oracle=lambda p,w:((p[0],p[0]),(F(2,5)-p[0],F(2,5)-p[0]))
        r=construct(oracle,[1,1],1,F(1,20),48)
        receipt=r['receipt'];self.assertEqual(receipt['status'],'conditional_certificate')
        self.assertLessEqual(receipt['lower'],F(2,5));self.assertGreaterEqual(receipt['upper'],F(2,5))
        self.assertGreaterEqual(receipt['policy'][0],F(2,5))
        self.assertEqual(receipt,certify(r['paths'],r['samples'],[1,1],1,F(1,20)))

    def test_two_step_conservation_and_heat_models(self):
        # Different physical equations check policy semantics and global bounds.
        def reaction(p,w):
            x=F(0)
            for u in p:x=x+u*(1-x)/2-x/4
            j=sum(p);g=F(9,20)-x
            return ((j,j),(g,g))
        def heat(p,w):
            hot=cold=F(0)
            for u in p:hot,cold=hot/2+cold/4+u/2,hot/4+cold/2
            j=sum(p);g=F(1,10)-cold;limit=hot-F(3,10)
            return ((j,j),(g,g),(limit,limit))
        for oracle,L,optimum in [(reaction,[2,F(7,8)],F(9,10)),(heat,[2,F(1,8),F(3,4)],F(4,5))]:
            r=construct(oracle,L,2,F(1,5),72);receipt=r['receipt']
            self.assertEqual(receipt['status'],'conditional_certificate')
            self.assertLessEqual(receipt['lower'],optimum)
            self.assertGreaterEqual(receipt['upper'],optimum)
            exact=oracle(receipt['policy'],F(0))
            self.assertTrue(all(hi<=0 for _,hi in exact[1:]))
            self.assertEqual(receipt,certify(r['paths'],r['samples'],L,2,F(1,5)))

    def test_bounded_exact_inputs(self):
        with self.assertRaises(ValueError):box_for('x',2)
        with self.assertRaises(ValueError):box_for('',9)
        with self.assertRaises(ValueError):certify([''],[Sample((.5,),((0,0),))],[1],1,1)
        with self.assertRaises(ValueError):certify([''],[Sample((F(0),),((0,0),))],[-1],1,1)


if __name__=='__main__':unittest.main()
