import unittest
import numpy as np
from experiments.v6_headroom import Trace,Scaling,Model,compile_state,forecast,features,physical_system,fit,blocks


class HeadroomTests(unittest.TestCase):
    def test_independent_recurrence_and_companion_match_all_orders(self):
        rng=np.random.default_rng(13)
        for p in (1,2,3,4,6,8):
            for k in (1,2):
                norm=Scaling(20.,3.,np.zeros(k),np.ones(k))
                coef=rng.normal(0,.03,p*(1+k)+1);coef[0]=.5
                m=compile_state(Model(p,norm,coef=coef))
                y=rng.normal(20,3,12);past=rng.normal(size=(12,k));future=rng.normal(size=(32,k))
                ordinary=forecast(m,y,past,future)
                compiled=forecast(m,y,past,future,True)
                np.testing.assert_allclose(ordinary,compiled,atol=1e-12,rtol=1e-12)

    def test_no_future_output_argument_and_no_history_aliasing(self):
        norm=Scaling(0.,1.,np.zeros(2),np.ones(2))
        m=compile_state(Model(2,norm,coef=np.array([.5,.1,.2,.1,.05,-.03,0.])))
        y=np.arange(6,dtype=float);u=np.ones((5,2));future=np.zeros((16,2))
        ycopy=y.copy();ucopy=u.copy()
        pred=forecast(m,y,u,future)
        np.testing.assert_array_equal(y,ycopy);np.testing.assert_array_equal(u,ucopy)
        self.assertEqual(pred.shape,(16,))
        with self.assertRaises(ValueError):forecast(m,np.array([1.]),u,future)

    def test_latent_linear_system_has_finite_history_closure(self):
        rng=np.random.default_rng(99);A=np.diag((.1,.4,.7))
        B=np.array([[1.,.2],[.1,.8],[.4,-.3]]);C=np.array([1.,.4,-.2])
        u=rng.normal(size=(600,2));x=np.zeros(3);ys=[C@x]
        for current in u:
            x=A@x+B@current;ys.append(C@x)
        trace=Trace(u,np.asarray(ys));norm=Scaling(0.,1.,np.zeros(2),np.ones(2))
        X,target=features(trace,3,norm)
        coef=np.linalg.lstsq(X[:400],target[:400],rcond=None)[0]
        self.assertLess(float(np.max(np.abs(X[400:]@coef-target[400:]))),1e-12)
        m=compile_state(Model(3,norm,coef=coef))
        prediction=forecast(m,trace.y[:451],trace.u[:450],trace.u[450:466],True)
        np.testing.assert_allclose(prediction,trace.y[451:467],atol=1e-12)

    def test_hidden_physics_is_confined_to_generator(self):
        train,select,cal,test,meta=physical_system(62001)
        self.assertEqual(meta['trajectory_count'],10)
        self.assertEqual(meta['simulator_transitions'],1920)
        self.assertGreaterEqual(meta['dimension'],3)
        self.assertTrue(all(v>0 for v in meta['generator_eigenvalues']))
        ordinary=fit(train,select,'stateful_arx')
        privileged=fit(train,select,'privileged_history_state')
        self.assertEqual(ordinary.p,privileged.p)
        np.testing.assert_array_equal(ordinary.coef,privileged.coef)
        np.testing.assert_allclose(blocks(ordinary,test),blocks(privileged,test,True),atol=1e-11)
        self.assertEqual(ordinary.work['candidate_fits'],24)
        self.assertEqual(ordinary.work['normal_equation_work_proxy'],privileged.work['normal_equation_work_proxy'])


if __name__=='__main__':unittest.main()
