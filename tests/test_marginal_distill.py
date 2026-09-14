from fractions import Fraction as F
import unittest

from experiments.marginal_distill import candidate_seeds
from experiments.marginal_orbit_quartic import symmetry_group
from experiments.marginal_symbolic import encode


class DistillTests(unittest.TestCase):
    def test_factor_and_orbit_sources_give_same_normalized_direction(self):
        word=((0,0),)
        factor={'modes':4,'operator_degree':4,'denominator':10,
                'blocks':[{'words':[word],'factor':[[2]]}]}
        orbit={'modes':4,'operator_degree':4,'permutations':symmetry_group(4),
               'orbit_squares':[{'polynomial':encode({word:F(1,5)}),'weight':'1'}]}
        self.assertEqual(candidate_seeds(factor),[{word:F(1)}])
        self.assertEqual(candidate_seeds(factor),candidate_seeds(orbit))
        orbit['permutations']=[[0,1,2,3]]
        with self.assertRaises(ValueError):candidate_seeds(orbit)


if __name__=='__main__':unittest.main()
