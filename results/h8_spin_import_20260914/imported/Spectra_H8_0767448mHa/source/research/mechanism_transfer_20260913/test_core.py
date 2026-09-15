"""Independent exact sixteen-mode control and preserved refusal/binding gates."""
from copy import deepcopy
from fractions import Fraction as F
from pathlib import Path
import json
import unittest

from experiments.marginal_symbolic import encode, mono
from research.molecular_collective_20260913.core import digest
from research.mechanism_transfer_20260913.core import replay


def control():
    h={((1,i),(0,i)):F(-1 if i<8 else 1) for i in range(16)}
    data={'modes':16,'particles':8,'hamiltonian':encode(h)}
    tail={'kind':'molecular_density_tail_v1','fixture_sha256':digest(data),
        'center_number':True,'factors':[],'alpha':'0','beta':'0'}
    blocks=[{'name':str(i),'words':[((1 if i<8 else 0,i),)],'factor':[[1]]} for i in range(16)]
    cert={'kind':'mechanism_spin_subspace_v1','fixture_sha256':digest(data),'tail_sha256':digest(tail),
        'denominator':1,'b':'-9','number_multiplier':[],'base_blocks':blocks,
        'anti_blocks':[{'name':'linear','generators':[[-1,-1,-1,0]],'directions':[[2]],
            'direction_denominator':2,'factor':[[1]]}]}
    return data,tail,cert


class TransferChecker(unittest.TestCase):
    def test_sixteen_mode_control_has_known_bound_and_zero_residual(self):
        d,t,c=control(); r=replay(d,t,c)
        self.assertEqual(F(r['original_lower_Ha']),-9)
        self.assertEqual(F(r['retained']['residual_l1']),0)
        # Eight negative one-particle levels give exact E0=-8; the added
        # {a0^dagger,a0}=I intentionally weakens the lower by exactly one.
        self.assertLessEqual(F(r['original_lower_Ha']),sum([-1]*8))

    def test_wrong_inputs_and_malformed_integer_factors_are_refused(self):
        d,t,c=control()
        for edit in ('fixture','direction','denominator','mode'):
            cc=deepcopy(c)
            if edit=='fixture': cc['fixture_sha256']='bad'
            if edit=='direction': cc['anti_blocks'][0]['directions']=[[True]]
            if edit=='denominator': cc['anti_blocks'][0]['direction_denominator']=0
            if edit=='mode': cc['anti_blocks'][0]['generators'][0][-1]=16
            with self.assertRaises(ValueError): replay(d,t,cc)
        dd={**d,'modes':18}; tt={**t,'fixture_sha256':digest(dd)}
        cc={**c,'fixture_sha256':digest(dd),'tail_sha256':digest(tt)}
        with self.assertRaises(ValueError): replay(dd,tt,cc)

    def test_inherited_h6_lower_is_identical(self):
        root=Path(__file__).resolve().parents[2]
        base=root/'results/molecular_collective_20260913/campaign/h6'
        source=root/'results/spin_completion_20260913/campaign/adaptive/round_4_batch_2'
        d=json.loads((base/'fixture.json').read_text()); t=json.loads((base/'rank_10/tail.json').read_text())
        c=json.loads((source/'certificate.json').read_text()); c['kind']='mechanism_spin_subspace_v1'
        r=replay(d,t,c); old=json.loads((source/'receipt.json').read_text())['accepted']
        self.assertEqual(F(r['original_lower_Ha']),F(old['original_lower_Ha']))


if __name__=='__main__': unittest.main()
