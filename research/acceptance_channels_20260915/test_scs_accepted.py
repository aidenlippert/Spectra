import unittest
from types import SimpleNamespace
import numpy as np
from scipy import sparse
import scs
from research.acceptance_channels_20260915.scs_accepted import accepted_problem
from research.acceptance_channels_20260915.scs_control import unpack


class AcceptedSCSTest(unittest.TestCase):
    def test_full_L1_conic_problem_matches_known_restricted_bound(self):
        # Only diagonal squares: the off-diagonal residual cannot be removed.
        # Max b-|2-b-q0|-|1-b-q1|-2|.2| = 0.6.
        op=SimpleNamespace(free=sparse.csc_matrix([[1.],[0.],[1.]]),
            M=[sparse.csr_matrix([[1.],[0.],[0.]]),sparse.csr_matrix([[0.],[0.],[1.]])],
            Q=[np.zeros((1,1)),np.zeros((1,1))],rhs=np.array([2.,.2,1.]))
        K=sparse.eye(3,format='csr');weights=np.array([1.,2.,1.])
        data,cone=accepted_problem(op,np.array([0]),K,weights)
        result=scs.SCS(data,cone,eps_abs=1e-10,eps_rel=1e-10,verbose=False).solve()
        b=-result['y'][0];q=result['y'][10:12]
        residual=op.rhs-np.array([q[0]+b,0,q[1]+b])
        self.assertAlmostEqual(b-weights@abs(residual),.6,places=7)
        self.assertAlmostEqual(result['info']['pobj'],.6,places=7)
        self.assertGreaterEqual(q.min(),-1e-8)

    def test_full_Gram_returns_known_small_matrix_eigenvalue(self):
        H=np.array([[2.,.2],[.2,1.]])
        op=SimpleNamespace(free=sparse.csc_matrix([[1.],[0.],[1.]]),
            M=[sparse.csr_matrix([[1.,0.,0.,0.],[0.,.5,.5,0.],[0.,0.,0.,1.]])],
            Q=[np.zeros((2,2))],rhs=np.array([2.,.2,1.]))
        data,cone=accepted_problem(op,np.array([0]),sparse.eye(3,format='csr'),np.array([1.,2.,1.]))
        result=scs.SCS(data,cone,eps_abs=1e-10,eps_rel=1e-10,verbose=False).solve()
        b=-result['y'][0];Q=unpack(result['y'][10:],2)
        residual=H-Q-b*np.eye(2)
        lower=b-abs(residual).sum()
        self.assertAlmostEqual(lower,np.linalg.eigvalsh(H)[0],places=7)


if __name__=='__main__':unittest.main()
