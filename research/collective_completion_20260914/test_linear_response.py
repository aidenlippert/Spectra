from fractions import Fraction as F
from itertools import combinations
import unittest
import numpy as np
from experiments.marginal_symbolic import mono,add,scale,product,adj,canonical,word_product
from experiments.marginal_coefficient import gram_map,dagger
from research.collective_completion_20260914.linear_response import embed,embed_paired_seed

class LinearResponseTests(unittest.TestCase):
    def test_cubic_cancellation_and_full_cross_terms(self):
        cubic=[((1,0),(0,2),(0,1)),((1,1),(0,3),(0,0))]
        linear=[((0,0),),((0,2),)];minus=cubic+linear;plus=[dagger(w) for w in minus]
        slots0=[0,1,2,3];slots1=[0,1,4,5];size=6
        rows=[tuple((1,i) for i in a)+tuple((0,i) for i in b) for k in range(4)
              for a in combinations(range(4),k) for b in combinations(range(4),k) if a<=b]
        lookup={w:i for i,w in enumerate(rows)}
        M=embed(gram_map(minus,lookup),slots0,size)+embed(gram_map(plus,lookup),slots1,size,True)
        self.assertEqual(M[[i for i,w in enumerate(rows) if len(w)==6],:].nnz,0)
        R=np.array([[2,-3,1,4,-2,3],[1,0,3,2,5,-1]],dtype=np.int64);Q=R.T@R
        symbolic={}
        for row in R:
            for words,slots in [(minus,slots0),(plus,slots1)]:
                p=add(*(mono(w,F(int(row[s]))) for w,s in zip(words,slots)))
                symbolic=add(symbolic,product(canonical(adj(p)),p))
        self.assertLessEqual(max(map(len,symbolic)),4)
        result=np.asarray(M@Q.ravel())
        self.assertEqual(result.tolist(),[float(symbolic.get(w,0)) for w in rows])
        # Independent linear coefficients alter quartic coefficients. Tying
        # the two responses would remove these permitted reconstruction terms.
        R2=R.copy();R2[:,slots1[2:]]=R2[:,slots0[2:]]
        d=np.asarray(M@((R.T@R-R2.T@R2).ravel()))
        self.assertGreater(sum(abs(d[i]) for i,w in enumerate(rows) if len(w)==4),0)
    def test_adjoint_cubic_identity(self):
        words=[((1,0),(0,2),(0,1)),((1,2),(0,3),(0,0))]
        for a in words:
            for b in words:
                p=add(dict(word_product(dagger(a),b)),dict(word_product(b,dagger(a))))
                self.assertLessEqual(max(map(len,p),default=0),4)
    def test_existing_paired_factors_survive_response_embedding(self):
        cubic=[((1,0),(0,2),(0,1)),((1,1),(0,3),(0,0))]
        words=cubic+[((0,0),),((0,2),)]
        V=np.array([[1,2],[-3,4],[5,-6],[7,8]],dtype=float)
        pair={'joint_size':6,'cubic_size':2,'slot_indices':[[0,1,2,3],[0,1,4,5]]}
        J=embed_paired_seed(words,V,words,pair)
        for old,new in zip(V.T,J[:,:2].T):
            for ws,slots in [(words,pair['slot_indices'][0]),([dagger(w) for w in words],pair['slot_indices'][1])]:
                before=add(*(mono(w,F(int(v))) for w,v in zip(ws,old)))
                after=add(*(mono(w,F(int(new[s]))) for w,s in zip(ws,slots)))
                self.assertEqual(canonical(before),canonical(after))
        np.testing.assert_array_equal(J[2:,2:],np.eye(4))
if __name__=='__main__':unittest.main()
