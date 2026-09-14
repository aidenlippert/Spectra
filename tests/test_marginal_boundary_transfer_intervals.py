import unittest
from fractions import Fraction as F

from experiments.marginal_boundary_transfer import (
    compile_state, contract, enclose, _dyadic, _interval_product,
)
from tests.test_marginal_boundary_transfer import act
from experiments.marginal_local_hubbard_block import _terms


VALENCE = {0x99: 1, 0x96: -1, 0x69: -1, 0x66: 1}


def block_state(q):
    state = {0: 1}
    for block in range(q):
        state = {s | (v << (8 * block)): a * b
                 for s, a in state.items() for v, b in VALENCE.items()}
    return state


class BoundaryTransferIntervalTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        moved=act(_terms(4,4,1),VALENCE)
        seed={s:8*VALENCE.get(s,0)-moved.get(s,0) for s in VALENCE.keys()|moved.keys()}
        cls.compiled = compile_state({s:a for s,a in seed.items() if a}, 4)

    def test_small_exact_contracts_are_enclosed(self):
        for a, b in ((F(1, 5), 0), (F(1, 5), F(1, 7))):
            for blocks in (2, 3, 7):
                with self.subTest(a=a, b=b, blocks=blocks):
                    exact = F(contract(self.compiled, a, b, blocks)['energy'])
                    result = enclose(self.compiled, a, b, blocks, 96)
                    self.assertTrue(result['accepted'])
                    self.assertLessEqual(F(result['energy_lower']), exact)
                    self.assertLessEqual(exact, F(result['energy_upper']))
                    self.assertLessEqual(F(result['width_per_site']), F(1, 10**20))

    def test_large_linear_valence_formula(self):
        # For the linear filter, the direct block energy is zero and each
        # merge contributes -3/13 at eta=1/5.
        result = enclose(compile_state(VALENCE,4), F(1, 5), 0, 125000, 96)
        expected = -F(3, 13) * 124999
        self.assertTrue(result['accepted'])
        self.assertLessEqual(F(result['energy_lower']), expected)
        self.assertLessEqual(expected, F(result['energy_upper']))
        self.assertEqual(result['sites'], 500000)

    def test_signed_matrix_interval_product_and_scale(self):
        bits = 96
        for scale in (F(1,2**60),1,2**120):
            a=[[F(1,3)*scale,F(-2,5)*scale],[F(-7,13)*scale,F(-5,7)*scale]]
            b=[[F(-7,11),F(2,3)],[F(3,13),F(-2,5)]]
            product=_interval_product(_dyadic(a,bits),_dyadic(b,bits),bits)
            unit=F(2)**(product[2]-bits)
            for i in range(2):
                for j in range(2):
                    lo,hi=product[0][i][j]*unit,product[1][i][j]*unit
                    exact=sum(a[i][k]*b[k][j] for k in range(2))
                    self.assertLessEqual(lo,exact)
                    self.assertLessEqual(exact,hi)

    def test_zero_norm_filter_is_refused(self):
        # Contact occupations are up/empty. h^2 acts as identity, so
        # I-h^2 annihilates this product at its only cut.
        compiled=compile_state({108:1},4)
        with self.assertRaises(ValueError):contract(compiled,0,-1,2)
        with self.assertRaises(ValueError):enclose(compiled,0,-1,2,96)

    def test_precision_and_size_refusals(self):
        with self.assertRaises(ValueError): enclose(self.compiled, 0, 0, 2, 63)
        with self.assertRaises(ValueError): enclose(self.compiled, 0, 0, 2, 513)
        with self.assertRaises(ValueError): enclose(self.compiled, 0, 0, 0, 96)
        with self.assertRaises(ValueError): enclose(self.compiled, 0, 0, 125000001, 96)
        with self.assertRaises(ValueError): enclose(self.compiled, 0.1, 0, 2, 96)


if __name__ == '__main__':
    unittest.main()
