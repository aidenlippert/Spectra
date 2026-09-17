"""Numerical proposer tests; run separately from accepting-checker tests."""
from fractions import Fraction as F
import unittest
import numpy as np
from research.intervention_reduction_20260916.predict import stable_step,schedule_from_trajectory,chebyshev_power_coefficients


class StableProposalTests(unittest.TestCase):
    def test_high_frequency_propagation_against_analytic_solution(self):
        a=np.array([[0.,0.],[0.,100.]])
        step=stable_step(a,F(1,4))
        self.assertLessEqual(float(step)*100,2)
        state=np.array([1.,1.],dtype=complex)/np.sqrt(2)
        remaining=F(1)
        while remaining:
            h=min(step,remaining);remaining-=h
            term=state.copy();candidate=term.copy()
            for k in range(1,19):
                term=(-1j*float(h)/k)*(a@term);candidate+=term
            state=candidate
        expected=np.array([1.,np.exp(-100j)])/np.sqrt(2)
        self.assertLess(np.linalg.norm(state-expected),1e-8)

    def test_zero_generator_keeps_requested_step(self):
        self.assertEqual(stable_step(np.zeros((2,2)),F(1,4)),F(1,4))

    def test_schedule_preserves_switch_time_and_controls(self):
        t={'segments':[{'duration_atomic_time':h,'controls_Ha':[u]}
           for h,u in [('1/4','1/100'),('3/4','1/100'),('1/2','-1/100')]]}
        self.assertEqual(schedule_from_trajectory(t),[(F(1),[F(1,100)]),(F(1,2),[F(-1,100)])])

    def test_integer_interpolation_against_analytic_high_frequency(self):
        count=25;nodes=(1+np.cos((np.arange(count)+.5)*np.pi/count))/2
        values=np.exp(-10j*nodes)[:,None];den=10**14
        power=chebyshev_power_coefficients(values,den)
        # Rational evaluation prevents large power coefficients from hiding error
        # in cancellation at theta near one.
        for theta in (F(0),F(1,7),F(1,2),F(9,10),F(1)):
            z=complex(sum(F(row[0][0],den)*theta**i for i,row in enumerate(power)),
                      sum(F(row[0][1],den)*theta**i for i,row in enumerate(power)))
            self.assertLess(abs(z-np.exp(-10j*float(theta))),1e-10)


if __name__=='__main__':
    unittest.main()
