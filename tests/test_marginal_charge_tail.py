import copy
from fractions import Fraction as F
import json
from pathlib import Path
import unittest
from unittest.mock import patch
from experiments.marginal_charge_tail import envelope,replay
from experiments.marginal_coherent_tree import CoherentCharge
from tests.test_marginal_charge_spin import fixed_sign_fixture
from tests.test_marginal_coherent_tree import direct_rows


class ChargeTailTests(unittest.TestCase):
    def test_envelope_is_below_every_physical_tail_row(self):
        data=fixed_sign_fixture();data['minimum_doublons']=2
        physical=direct_rows(CoherentCharge(data))
        with patch('experiments.marginal_determinant_tree.DeterminantOracle.action',side_effect=AssertionError('No action')):
            o,poly,_,_,cost=envelope(data)
        count=0
        for state,value in physical.items():
            if sum(((state>>(2*i))&3)==3 for i in range(o.sites))<2:continue
            lower=sum((v for mask,v in poly.items() if state&mask==mask),F(0))
            self.assertLessEqual(lower,value);count+=1
        self.assertEqual(count,6);self.assertEqual(cost['determinant_actions'],0)

    def test_exact_bracket_and_refusals(self):
        data=fixed_sign_fixture();data.update(kind='weighted_charge_tail_bracket_v1',minimum_doublons=2,
            witness_state=15,b='0',positive_indicators=[],charge_indicators=[],number_multipliers=[[],[]])
        r=replay(data);self.assertLessEqual(F(r['lower']),F(r['upper']))
        for mutate in [lambda c:c.update(witness_state=0),lambda c:c.update(minimum_doublons=0),
                       lambda c:c.update(positive_indicators=[{'required':0,'occupied':0,'weight':'-1'}]),
                       lambda c:c.update(number_multipliers=[[{'mask':15,'coefficient':'1'}],[]])]:
            bad=copy.deepcopy(data);mutate(bad)
            with self.assertRaises(ValueError):replay(bad)

    def test_saved_h6_bracket_refutes_requested_envelope_threshold(self):
        root=Path(__file__).resolve().parents[1]/'results/marginal_h6/charge_tail/bracket'
        c=json.loads((root/'certificate.json').read_text())
        with patch('experiments.marginal_determinant_tree.DeterminantOracle.action',side_effect=AssertionError('No action')):
            r=replay(c)
        self.assertEqual(r,json.loads((root/'receipt.json').read_text()))
        self.assertLess(F(r['upper']),F('-6.264'));self.assertLess(F(r['width']),F(1,10**9))
        self.assertEqual(r['explicit_upper_witnesses'],1)


if __name__=='__main__':unittest.main()
