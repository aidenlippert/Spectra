import unittest
from fractions import Fraction as F
from itertools import product
from decimal import Decimal,localcontext
from .collective_control import real_pauli,comm,controls,trig_interval,pi_interval,construct,replay
from .test_tensor_operator import literal_action
from experiments.marginal_hunt_car import scale
from .collective_obstruction import selected_action,all_angle_floor
from .transfer_collective import choose_pair
from research.molecular_collective_20260913.core import digest
import copy

class TestCollectiveControl(unittest.TestCase):
 def test_exact_control_algebra(self):
  d,w,a=controls(0,2)
  self.assertEqual(comm(w,d),scale(a,2));self.assertEqual(comm(w,a),scale(d,2))
  for pair in [(1,2),(0,0),(-2,2),(0,3)]:
   with self.assertRaises(ValueError):controls(*pair)
 def test_selected_action_against_independent_oracle(self):
  words=[(),((1,0),(0,3)),((1,2),(1,0),(0,3),(0,1)),((0,2),(1,2)),((1,0),(1,0))]
  for word in words:
   for ket in range(16):
    out,sign=literal_action(word,ket)
    self.assertEqual(selected_action({word:F(7,13)},ket),{out:F(7*sign,13)} if sign else {})
 def test_generalized_certificate_and_refusals(self):
  data={'modes':4,'particles':2,'hamiltonian':[]}
  cert={'kind':'integer_charge_mps_v1','fixture_sha256':digest(data),'modes':4,'particles':2,
        'spin_counts':[1,1],'denominator':1,'bond_charges':[[[0,0]],[[1,0]],[[1,1]],[[1,1]],[[1,1]]],
        'tensors':[[[0,s,0,1]] for s in (1,1,0,0)]}
  pair,occ,_=choose_pair(data,cert)
  self.assertEqual(pair,(0,2));self.assertEqual(occ,[F(2),F(0)])
  a=construct(data,cert,pair=pair,max_amplitude=F(1,2))
  self.assertTrue(a['accepted']);self.assertEqual(a['HD_norm_bound_Ha'],'0')
  self.assertTrue(replay(data,cert,a,F(1,2),expected_pair=pair)['accepted'])
  with self.assertRaises(ValueError):replay(data,cert,a,F(1))
  with self.assertRaises(ValueError):replay(data,cert,a,F(1,2),expected_pair=(2,0))
  with self.assertRaises(ValueError):construct(data,cert,pair=(0,4))
  bad=copy.deepcopy(a);bad['uncertainties']['per_control_pointwise_Ha']='0'
  with self.assertRaises(ValueError):replay(data,cert,bad,F(1,2),expected_pair=pair)
  bad=copy.deepcopy(cert);bad['tensors'][0][0][3]=1.0
  with self.assertRaises(ValueError):construct(data,bad,pair=pair)
 def test_all_angle_envelope_floor(self):
  # Exact conservative sample values below the H8 norm witnesses.
  floor,parts=all_angle_floor(F(139,100),F(61,100),F(76,100),F(1,2))
  self.assertEqual(floor,0);self.assertGreater(parts[1],1);self.assertGreater(parts[2],1)
  # Vanishing drift must never imply an obstruction to a pi pulse.
  floor,_=all_angle_floor(F(139,100),F(0),F(0),F(1,2))
  self.assertLessEqual(floor,F(-139,100))
 def test_pauli_conversion_against_independent_action(self):
  mats={'I':(1,0,0,1),'X':(0,1,1,0),'J':(0,-1,1,0),'Z':(1,0,0,-1)}
  words=[(),((1,0),(0,3)),((1,2),(1,0),(0,3),(0,1)),((0,2),(1,2)),((1,0),(1,0))]
  for word in words:
   terms=real_pauli({word:F(7,13)},4)
   for ket in range(16):
    out,sg=literal_action(word,ket)
    for bra in range(16):
     value=F()
     for w,c in terms.items():
      v=c
      for site,s in enumerate(w):v*=mats[s][2*((bra>>site)&1)+((ket>>site)&1)]
      value+=v
     self.assertEqual(value,F(7*sg,13) if bra==out else F())
 def test_trig_and_pi_independent_decimal_series(self):
  with localcontext() as ctx:
   ctx.prec=150;x=Decimal(25)/8
   for sine in [False,True]:
    t=x if sine else Decimal(1);s=t
    for k in range(1,140):
     a=2*k+int(sine);t*=-x*x/Decimal(a*(a-1));s+=t
    lo,hi=trig_interval(F(25,8),sine)
    self.assertLessEqual(Decimal(lo.numerator)/lo.denominator,s)
    self.assertGreaterEqual(Decimal(hi.numerator)/hi.denominator,s)
  lo,hi=pi_interval();self.assertTrue(lo>F(314159265358979323846,10**20));self.assertTrue(hi<F(314159265358979323847,10**20))

if __name__=='__main__':unittest.main()
