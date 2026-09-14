import unittest
from fractions import Fraction as F
from itertools import repeat
from experiments.v10_exp_readout import evaluate_modes,evaluate_mode_map,times,norm,plus,scale


def independent_series(z,n=100):
    # No scaling, squaring or grid rounding. This much longer series encloses
    # the analytic value independently of the readout implementation.
    term=(F(1),F(0));s=term
    for k in range(1,n+1):term=scale(times(term,z),F(1,k));s=plus(s,term)
    r=norm(z)/F(n+1)
    return s,norm(term)*r/(1-r)


class ReadoutTests(unittest.TestCase):
    def test_exact_and_cancellation(self):
        m,e,_=evaluate_modes([(0,0,0,3,-2)],0,F(1,1000))
        self.assertEqual(m,(3,-2));self.assertEqual(e,0)
        m,e,c=evaluate_modes([(0,9,1,10**40,0),(0,9,0,-10**40,0)],1,F(1,10**6))
        self.assertEqual(m,(0,0));self.assertEqual(e,0);self.assertEqual(c['unique_exponents'],0)

    def test_independent_complex_enclosures(self):
        for z in [(F(-1),F(0)),(F(-3),F(7)),(F(0),F(8)),(F(-1,10),F(2,3))]:
            m,e,c=evaluate_modes([(z[0],z[1],0,1,0)],1,F(1,10**8))
            other,other_error=independent_series(z)
            self.assertLessEqual(norm((m[0]-other[0],m[1]-other[1]))+other_error,e)
            self.assertLessEqual(e,F(1,10**8));self.assertGreater(c['squarings'],0)

    def test_shared_exponents_total_budget(self):
        rows={'P':[(0,j,0,j,1) for j in range(1,10)],
              'M':[(0,j,1,2,-j) for j in range(1,10)]}
        mids,e,c=evaluate_mode_map(rows,F(1,2),F(1,100000))
        self.assertEqual(c['unique_exponents'],9);self.assertEqual(c['input_terms'],18)
        self.assertLessEqual(e,F(1,100000));self.assertEqual(set(mids),{'P','M'})

    def test_refusals_are_bounded(self):
        for rows,tol in [([(1,0,0,1,0)],F(1,10)), ([(0,0,True,1,0)],F(1,10)),
                         ([(0,0,0,.5,0)],F(1,10)), ([(0,0,33,1,0)],F(1,10)),
                         ([(0,0,0,1,0)],0), ([(0,1<<30,0,1,0)],F(1,10))]:
            with self.assertRaises(ValueError):evaluate_modes(rows,1,tol)
        with self.assertRaises(ValueError):evaluate_modes(repeat((0,0,0,1,0)),1,F(1,10))
        with self.assertRaises(ValueError):evaluate_modes([],F(1,1<<8193),F(1,10))


if __name__=='__main__':unittest.main()
