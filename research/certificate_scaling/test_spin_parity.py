import unittest
from fractions import Fraction
from research.certificate_scaling.spin_parity import spin_parity

class SpinParityTests(unittest.TestCase):
    def test_span_bruteforce(self):
        for m in (4,6):
            p={((1,0),(0,2)):1, ((1,1),(0,3)):2}
            b=spin_parity(p,m)
            got={0}
            for x in b: got |= {y^x for y in tuple(got)}
            expected={x for x in range(1<<(m//2)) if all((x & r).bit_count()%2==0 for r in (3,))}
            lifted={sum(((x>>s)&1)*((1<<(2*s))|(1<<(2*s+1))) for s in range(m//2)) for x in expected}
            self.assertEqual(got,lifted)

    def test_intersection_combination(self):
        # Constraints x0=x1 and x1=x2; neither individual generator survives filtering.
        p={((1,0),(0,2)):1, ((1,2),(0,4)):1}
        self.assertEqual(len(spin_parity(p,6)),1)

    def test_density_quartic_and_zero_coefficients(self):
        p={((1,0),(0,0)):Fraction(2),((1,0),(1,2),(0,0),(0,2)):Fraction(3),((1,0),(0,2)):Fraction(0)}
        self.assertEqual(spin_parity(p,4),(3,12))

    def test_invalid(self):
        for p,m in [({((1,0),(0,1)):1.0},4),({((0,1),(1,0)):1},4),({},3)]:
            with self.assertRaises(ValueError): spin_parity(p,m)

if __name__=='__main__': unittest.main()
