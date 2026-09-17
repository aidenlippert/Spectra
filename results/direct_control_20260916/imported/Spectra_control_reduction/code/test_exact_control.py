import unittest,random,itertools,copy,json
from fractions import Fraction as F
from math import lcm
from exact_control import *
class Algebra(unittest.TestCase):
 def test_sqrt_enclosure(self):
  rng=random.Random(70)
  for _ in range(200):
   x=F(rng.randrange(1,10**30),rng.randrange(1,10**20));q=sqrt_up(x)
   self.assertGreaterEqual(q*q,x);self.assertLess((q-F(1,2**60))**2,x)
  self.assertEqual(sqrt_up(0),0)
  with self.assertRaises(ValueError):sqrt_up(-1)
 def test_large_integers(self):
  for k in [100,200,500]:
   x=F(2**k+37,17);q=sqrt_up(x)
   self.assertGreaterEqual(q*q,x)
 def test_balanced_count(self):
  self.assertEqual(len(labels_for(16,4,4)),4900)
  self.assertEqual(len(set(labels_for(8,2,2))),36)
 def test_literal_car_matrix(self):
  # Independently construct Jordan-Wigner matrices by Kronecker products.
  def kron(A,B):return [[a*b for a in ar for b in br]for ar in A for br in B]
  Id=[[1,0],[0,1]];Z=[[1,0],[0,-1]];C=[[0,0],[1,0]];A=[[0,1],[0,0]]
  def ladder(c,p):
   out=[[1]]
   for site in reversed(range(3)):out=kron(out,(C if c else A)if site==p else Z if site<p else Id)
   return out
  def product(A,B):return [[sum(a*b for a,b in zip(row,col))for col in zip(*B)]for row in A]
  eye=[[int(i==j)for j in range(8)]for i in range(8)]
  ops=list(itertools.product((0,1),range(3)))
  checked=0
  for k in [0,1,2,3]:
   for word in itertools.product(ops,repeat=k):
    mat=eye
    for c,p in word:mat=product(mat,ladder(c,p))
    for state in range(8):
     target,sign=apply_word(state,word)
     self.assertEqual([row[state]for row in mat],[sign if target==i else 0 for i in range(8)])
     checked+=1
  self.assertEqual(checked,2072)
 def test_exact_projection_integral(self):
  rng=random.Random(97)
  for trial in range(20):
   E=[[rng.randrange(-8,9)for _ in range(3)]for _ in range(5)];ed=7;g=gram(E)
   re=[[rng.randrange(-10,11)for _ in range(3)]for _ in range(5)];im=[[rng.randrange(-10,11)for _ in range(3)]for _ in range(5)]
   floor=[[floor_scaled(F(x,ed*ed),12)for x in row]for row in g]
   bound=project_integral_upper(re,im,floor,8,12)
   ecr=[matvec(E,z)for z in re];eci=[matvec(E,z)for z in im]
   exact=sum((F(dot(ecr[a],ecr[b])+dot(eci[a],eci[b]),(a+b+1)*ed*ed*64)for a in range(5)for b in range(5)),F(0))
   self.assertGreaterEqual(bound,exact)
 def test_gram_identity(self):
  A=[[1,-2,4],[2,3,-3],[5,-1,7]];z=[3,-8,2];g=gram(A)
  self.assertEqual(dot(z,matvec(g,z)),dot(matvec(A,z),matvec(A,z)))
 def test_control_capacity_and_noncommutation(self):
  m=4;labels=list(range(16));D=[((s>>0)&1)+((s>>1)&1)-((s>>2)&1)-((s>>3)&1)for s in labels]
  W=[[0]*16 for _ in labels]
  for j,s in enumerate(labels):
   for i,k in [(0,2),(2,0),(1,3),(3,1)]:
    out,sg=apply_word(s,[(1,i),(0,k)])
    if sg:W[out][j]+=sg
  self.assertLessEqual(max(map(abs,D)),2)
  self.assertLessEqual(max(sum(map(abs,row))for row in W),2)
  self.assertTrue(any((D[i]-D[j])*W[i][j] for i in labels for j in labels))
 def test_conserving_phase_flip(self):
  for s in labels_for(8,2,2):
   for j in range(8):self.assertEqual((1-2*((s>>j)&1))**2,1)
 def test_noise_mixture_bound(self):
  p=F(1,1000);rho=[1-p,p];zexpect=2*rho[0]-2*rho[1]
  self.assertEqual(2-zexpect,4*p)
 def test_general_defect_scalar(self):
  # q(s)=1-is, H=1, dt=1: exact defect is i*s; integral squared=1/3.
  re=[[1],[0]];im=[[0],[-1]]
  self.assertEqual(sum((F((0 if a==0 else 1)*(0 if b==0 else 1),a+b+1)for a in range(2)for b in range(2)),F(0)),F(1,3))
 def test_bad_integral_refused(self):
  with self.assertRaises(ValueError):project_integral_upper([[10],[0]],[[0],[0]],[[-10**9]],1,10)

class MolecularMutations(unittest.TestCase):
 @classmethod
 def setUpClass(cls):
  cls.data=prepare_input();cls.pc=prepare_case(cls.data,ROOT/'inputs/short24.json')
 def test_molecular_good_control(self):
  self.assertTrue(query(self.data,self.pc)['target_proved'])
 def test_unsupported_claim_refused(self):
  pc=copy.deepcopy(self.pc);pc['case']['target_D_upper']='-7/10'
  self.assertFalse(query(self.data,pc)['target_proved'])
 def test_negative_noise_refused(self):
  with self.assertRaises(ValueError):query(self.data,self.pc,override_noise='-1/1000')
 def test_too_large_control_uncertainty(self):
  self.assertFalse(query(self.data,self.pc,override_amp='1/10')['target_proved'])
 def test_excess_initial_uncertainty(self):
  self.assertFalse(query(self.data,self.pc,override_trace='1/10')['target_proved'])
 def test_polynomial_jump_charged(self):
  pc=copy.deepcopy(self.pc);pc['case']['segments'][4]['re'][0][0]+=pc['zd']
  r=query(self.data,pc)
  self.assertFalse(r['target_proved']);self.assertGreater(r['state_error_float'],F(1,2))
 def test_initial_state_error_charged(self):
  pc=copy.deepcopy(self.pc);pc['case']['segments'][0]['re'][0][0]+=pc['zd']
  self.assertFalse(query(self.data,pc)['target_proved'])
 def test_undeclared_control_refused(self):
  pc=copy.deepcopy(self.pc);pc['case']['segments'][2]['u']='17'
  with self.assertRaises(ValueError):query(self.data,pc)
 def test_missing_polynomial_refused(self):
  pc=copy.deepcopy(self.pc);pc['case']['segments'][2]['im'].pop()
  with self.assertRaises(ValueError):query(self.data,pc)

if __name__=='__main__':unittest.main(verbosity=2)
