from fractions import Fraction as F
import unittest
from unittest.mock import patch
from types import SimpleNamespace
from experiments.marginal_charge_ratio_dp import ratio_upper
from experiments.marginal_coherent_tree import CoherentCharge
from experiments.marginal_spin_constructor import spin_states
from tests.test_marginal_coherent_tree import fixture


class ChargeRatioDPTests(unittest.TestCase):
    def test_all_transition_ratios_match_independent_occupation_maxima(self):
        o=CoherentCharge(fixture());states=spin_states(o.oracle);checked=0
        for c,a in o.groups:
            delta=tuple(((c>>(2*i))&3).bit_count()-((a>>(2*i))&3).bit_count() for i in range(o.sites))
            for minimum in (0,1,2):
                with patch('experiments.marginal_determinant_tree.DeterminantOracle.action',side_effect=AssertionError('No actions')):
                    value,cost=ratio_upper(o,delta,c|a,a,minimum)
                candidates=[]
                for s in states:
                    if s&(c|a)!=a:continue
                    q=[((s>>(2*i))&3).bit_count()-1 for i in range(o.sites)]
                    if q.count(1)<minimum:continue
                    ratio=F(1)
                    for label,factor in o.local_factors:
                        before=after=1
                        for i,power in label:before*=q[i]**power;after*=(q[i]+delta[i])**power
                        ratio*=factor**(after-before)
                    candidates.append(ratio)
                self.assertEqual(value,max(candidates) if candidates else None)
                self.assertEqual(cost['occupation_endpoint_evaluations'],0);checked+=1
        self.assertGreater(checked,10)

    def test_empty_domain_and_invalid_inputs(self):
        o=CoherentCharge(fixture())
        value,_=ratio_upper(o,(0,)*4,o.all_bits,0,1);self.assertIsNone(value)
        value,cost=ratio_upper(o,(0,)*4,0,0,1);self.assertEqual(value,1);self.assertEqual(cost['memory_sites'],0)
        for delta,mask,bits,d in [((0,)*3,0,0,1),((0,)*4,0,1,1),((0,)*4,0,0,3),((0,)*4,0,0,True)]:
            with self.assertRaises(ValueError):ratio_upper(o,delta,mask,bits,d)
        wide=SimpleNamespace(sites=6,all_bits=(1<<12)-1,target=3,
                             local_factors=[(((0,1),(5,1)),F(2))])
        with self.assertRaises(ValueError):ratio_upper(wide,(1,0,0,0,0,-1),0,0,1)


if __name__=='__main__':unittest.main()
