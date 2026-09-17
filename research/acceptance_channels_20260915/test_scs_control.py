import unittest
from types import SimpleNamespace
import numpy as np
from scipy import sparse
import scs
from research.acceptance_channels_20260915.scs_control import pack,unpack,standard_problem


class SCSControlTest(unittest.TestCase):
    def test_real_SCS_problem_recovers_known_matrix_lower_bound(self):
        H=np.array([[2.,.2],[.2,1.]])
        op=SimpleNamespace(free=sparse.csc_matrix([[1.],[0.],[1.]]),
            M=[sparse.csr_matrix([[1.,0.,0.,0.],[0.,.5,.5,0.],[0.,0.,0.,1.]])],
            Q=[np.zeros((2,2))],rhs=np.array([2.,.2,1.]))
        data,cone=standard_problem(op,np.array([0]))
        result=scs.SCS(data,cone,eps_abs=1e-10,eps_rel=1e-10,verbose=False).solve()
        b=-result['y'][0];Q=unpack(result['y'][1:],2)
        np.testing.assert_allclose(Q+b*np.eye(2),H,atol=1e-8)
        self.assertAlmostEqual(b,np.linalg.eigvalsh(H)[0],places=7)
        np.testing.assert_allclose(unpack(pack(H),2),H,atol=1e-15)


if __name__=='__main__':unittest.main()
