from fractions import Fraction as F
import json
import unittest
import numpy as np
import scipy.sparse as sp
from experiments.marginal_asymmetric_polish import solve_weights,native_lp
from experiments.marginal_asymmetric_columns import prepare_candidates,price,center_dual
from experiments.marginal_asymmetric_adapt import exchange_projection
from experiments.marginal_reynolds import invariant_rows,reynolds_matrix
from experiments.marginal_symbolic import product,adj
from experiments.marginal_orbit_certificate import average_canonical_polynomial


class AsymmetricColumnTests(unittest.TestCase):
    def test_native_small_entry_threshold_preserves_equation(self):
        matrix=sp.csr_matrix([[1.,1e-10]]);rhs=np.ones(1);bounds=[(None,None),(100,100)]
        coarse=native_lp(np.array([-1.,0.]),matrix,rhs,bounds,1e-9)
        fine=native_lp(np.array([-1.,0.]),matrix,rhs,bounds,1e-12)
        self.assertTrue(coarse.success and fine.success)
        self.assertGreater(abs((matrix@coarse.x-rhs)[0]),1e-9)
        self.assertLess(abs((matrix@fine.x-rhs)[0]),1e-14)

    def test_centered_dual_respects_normalization_and_energy_face(self):
        data={'unit':np.array([1.,0.]),'a':np.array([[0.],[1.]]),'rhs':np.array([1.,0.]),'weights':np.ones(2)}
        dual,status=center_dual(data,[np.array([1.,0.])],1.)
        self.assertTrue(status['accepted'])
        json.dumps(status)
        self.assertTrue(np.allclose(dual,[1.,0.],atol=1e-6,rtol=0))

    def test_highs_dual_sign_prices_an_improving_column(self):
        data={'projection':sp.eye(1),'basis':[],'unit':np.ones(1),'a':np.zeros((1,0)),
              'weights':np.ones(1),'rhs':np.ones(1),'source':'scalar dual-sign test','polish_build_seconds':0}
        result,_=solve_weights(data,[('test',{})],[np.ones(1)])
        self.assertTrue(result.success)
        moment=-result.eqlin.marginals
        self.assertAlmostEqual(float(data['unit']@moment),1)
        self.assertLess(float(-moment[0]),0)
        expanded,status=solve_weights(data,[('test',{}),('test',{})],[np.ones(1),-np.ones(1)],max_weight=2)
        self.assertGreater(-expanded.fun,-result.fun)
        self.assertAlmostEqual(-expanded.fun,3.)
        self.assertEqual(status['capped_square_weights'],1)

    def test_priced_columns_match_exchange_average_and_allow_asymmetry(self):
        modes=6;rows=invariant_rows(modes,4);projection,_=exchange_projection(rows,modes)
        data={'modes':modes,'rows':rows,'projection':projection}
        blocks,maps=prepare_candidates(data);dual=np.random.default_rng(194).normal(size=projection.shape[0])
        proposals,count=price(blocks,maps,dual,8)
        self.assertGreater(count,0);different=False;full=reynolds_matrix(rows,modes)[0]
        for proposal in proposals:
            p=proposal['polynomial'];square=product(adj(p),p)
            exact=average_canonical_polynomial(square,[list(range(modes)),[(i+3)%modes for i in range(modes)]])
            expected=projection@np.array([float(exact.get(w,0)) for w in rows])
            self.assertTrue(np.allclose(proposal['column'],expected,atol=1e-12,rtol=0))
            self.assertAlmostEqual(float(dual@proposal['column']),proposal['violation'],places=10)
            symmetric=projection@full@np.array([float(square.get(w,0)) for w in rows])
            different|=bool(np.max(abs(symmetric-proposal['column']))>1e-5)
        self.assertTrue(different)


if __name__=='__main__':unittest.main()
