import unittest,json,time
from fractions import Fraction
from pathlib import Path
from spin_irrep import *
from experiments.marginal_symbolic import canonical,add,scale,product,mono

class HighestWeightChecks(unittest.TestCase):
 @classmethod
 def setUpClass(cls):
  cls.op=load_operator();cls.groups=cls.op.meta['groups'];cls.R,cls.next=raising(cls.groups)
 def test_all_raising_words_against_literal_CAR(self):
  S=add(*(mono(((1,i),(0,i+1))) for i in range(0,16,2)));count=0
  for g in range(42,58):
   for j,w in enumerate(self.groups[g]['words']):
    B=canonical(mono(tuple(map(tuple,w))));literal=add(product(S,B),scale(product(B,S),-1));proposed={}
    if g in self.R:
     col=self.R[g].getcol(j).tocoo();h=self.next[g]
     for i,v in zip(col.row,col.data):proposed=add(proposed,scale(canonical(mono(tuple(map(tuple,self.groups[h]['words'][i])))),int(v)))
    self.assertEqual(literal,proposed,(g,j));count+=1
  self.assertEqual(count,3840)
 def test_exact_kernel_and_rank(self):
  for g,U in self.R.items():
   if U.shape!=(112,368):continue
   self.assertTrue(np.array_equal((U@U.T).toarray(),3*np.eye(112,dtype=int)))
   K=kernel(U);K6=np.rint(6*K).astype(np.int64)
   self.assertEqual(K6.shape,(368,256));self.assertTrue(np.all(U@K6==0))
   gram=K6.T@K6;self.assertTrue(np.all(gram==np.diag(np.diag(gram))));self.assertTrue(np.all(np.diag(gram)>0))
 def test_exact_spin_commutator(self):
  for g in range(42,58):
   n=len(self.groups[g]['words']);A=np.zeros((n,n),dtype=np.int64)
   if g in self.R:A-=(self.R[g].T@self.R[g]).toarray()
   for h,target in self.next.items():
    if target==g:A+=(self.R[h]@self.R[h].T).toarray()
   w=self.groups[g]['words'][0];twom=sum((1 if c else -1)*(1 if p%2==0 else -1) for c,p in w)
   self.assertTrue(np.array_equal(A,twom*np.eye(n,dtype=np.int64)))
 def test_random_full_gram_embedding(self):
  op=load_operator();rng=np.random.default_rng(6623)
  for i,m in enumerate(op.members):
   if m[0]>=42:op.V[i]=np.eye(len(op.V[i]));op.identity[i]=True;op.Q[i]=np.zeros((len(op.V[i]),)*2)
  pos={m[0]:i for i,m in enumerate(op.members)}
  for repeat in range(2):
   qs=[]
   for i,v in enumerate(op.V):
    v0=rng.normal(size=(v.shape[1],2));qs.append(v0@v0.T)
   specs=construct(op,qs);new=sum((op.M[pos[g]]@((v@q@v.T).ravel()) for g,v,q,_ in specs),np.zeros(len(op.rhs)));old=op.A(qs)
   self.assertLess(np.linalg.norm(new-old)/np.linalg.norm(old),2e-13)
 def test_reject_wrong_sector_and_hamiltonian(self):
  from research.collective_completion_20260914.spin_screen import check_sector
  fixture=json.loads(Path(self.op.meta['fixture']).read_text());c=json.loads((OUT/'spin_polish003/export/certificate.json').read_text())
  wrong=json.loads(json.dumps(c));wrong['particles']=7
  with self.assertRaises(ValueError):check_sector(fixture,wrong)
  wrong=json.loads(json.dumps(c));wrong['hamiltonian'][0]['coefficient']='123'
  with self.assertRaises((ValueError,KeyError)):check_sector(fixture,wrong)
 def test_full_accepted_rational_interval(self):
  r=json.loads((OUT/'spin_polish003/exact_interval.json').read_text());u=json.loads((OUT/'upper_replay.json').read_text())
  self.assertEqual(Fraction(r['upper_Ha']),Fraction(u['upper_Ha']))
  self.assertEqual(Fraction(r['width_Ha']),Fraction(r['upper_Ha'])-Fraction(r['lower']))
  self.assertLessEqual(Fraction(r['width_Ha']),Fraction(1,625));self.assertEqual(r['many_body_states_enumerated'],0)

if __name__=='__main__':unittest.main(verbosity=2)
