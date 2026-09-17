import unittest
import numpy as np
from .charge_mps import compress,linear_combination,apply_mpo,structural_error
from .test_mps_numeric import dense

def charged_mps(seed=2):
    rng=np.random.default_rng(seed);n=4;q=[]
    for i in range(n+1):
        q.append([(a,b) for a in range(2) for b in range(2)
                  if a<=sum(j%2==0 for j in range(i)) and b<=sum(j%2==1 for j in range(i))
                  and 1-a<=sum(j%2==0 for j in range(i,n)) and 1-b<=sum(j%2==1 for j in range(i,n))])
    out=[]
    for i in range(n):
        t=np.zeros((len(q[i]),2,len(q[i+1])))
        for l,ql in enumerate(q[i]):
            for s in (0,1):
                qr=(ql[0]+(s if i%2==0 else 0),ql[1]+(s if i%2 else 0))
                for r,x in enumerate(q[i+1]):
                    if x==qr:t[l,s,r]=rng.normal()
        out.append(t)
    return out,q

class ChargeMPSTests(unittest.TestCase):
    def test_exact_compression_at_high_bond(self):
        a,q=charged_mps();c,cq,d=compress(a,q,8)
        np.testing.assert_allclose(dense(c),dense(a),atol=1e-12)
        self.assertEqual(structural_error(c,cq),0)
        self.assertEqual(cq[0],[(0,0)]);self.assertEqual(cq[-1],[(1,1)])
        self.assertLess(d['discarded_squared_norm'],1e-20)

    def test_linear_combination_with_different_bonds(self):
        a,q=charged_mps(3);b,bq=charged_mps(5)
        b[0]=np.concatenate([b[0],np.zeros((1,2,1))],axis=2)
        b[1]=np.concatenate([b[1],np.zeros((1,2,b[1].shape[2]))],axis=0)
        bq[1].append((0,0))
        out,oq=linear_combination([(a,q),(b,bq)],[1.25,-.4])
        np.testing.assert_allclose(dense(out),1.25*dense(a)-.4*dense(b),atol=1e-12)
        self.assertEqual(structural_error(out,oq),0)

    def test_nontrivial_virtual_mpo_application(self):
        a,q=charged_mps(5);I=np.eye(2);Z=np.diag([1.,-1.]);CP=np.array([[0.,0.],[1.,0.]]);AN=CP.T
        w0=np.zeros((1,2,2,2));w1=np.zeros((2,2,2,2));w2=w1.copy();w3=np.zeros((2,1,2,2))
        w0[0,0]=I;w0[0,1]=CP;w1[0,0]=I;w1[1,1]=Z;w2[0,0]=I;w2[1,1]=AN;w3[:,0]=I
        oq=[[(0,0)],[(0,0),(1,0)],[(0,0),(1,0)],[(0,0),(0,0)],[(0,0)]]
        out,outq=apply_mpo(a,q,[w0,w1,w2,w3],oq)
        want=(np.eye(16)+np.kron(np.kron(np.kron(CP,Z),AN),I))@dense(a)
        np.testing.assert_allclose(dense(out),want,atol=1e-12)
        self.assertEqual(structural_error(out,outq),0)
        c,cq,_=compress(out,outq,8)
        np.testing.assert_allclose(dense(c),want,atol=1e-12)

if __name__=='__main__':unittest.main()
