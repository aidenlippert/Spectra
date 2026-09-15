"""Exact differential and refusal checks for compiled integer Gram products."""
import copy
import random
import unittest
from experiments.marginal_symbolic import expand_squares as original
from research.nvidia_followup_20260915.flint_squares import expand_squares as compiled


class GramTests(unittest.TestCase):
    def test_exact_small_large_and_empty_integer_factors(self):
        rng=random.Random(715)
        words=[[(0,0)],[(0,1)],[(1,2),(0,2),(0,0)],[(1,1),(0,1),(0,2)]]
        for rows in (0,1,5,12):
            for magnitude in (3,10**9,10**35):
                blocks=[{'words':words,'factor':[[rng.randrange(-magnitude,magnitude+1) for _ in words] for _ in range(rows)]}]
                self.assertEqual(original(blocks,10**8,4),compiled(blocks,10**8,4))

    def test_original_refusal_paths(self):
        base=[{'words':[[(0,0)],[(0,1)]],'factor':[[1,2]]}]
        malformed=[]
        a=copy.deepcopy(base);a[0]['factor'][0][0]=True;malformed.append(a)
        a=copy.deepcopy(base);a[0]['factor'][0].append(1);malformed.append(a)
        a=copy.deepcopy(base);a[0]['words'][0]=[(True,0)];malformed.append(a)
        a=copy.deepcopy(base);a[0]['words'][0]=[(1,0)];malformed.append(a)
        a=copy.deepcopy(base);a[0]['words'][0]=[(0,9)];malformed.append(a)
        malformed.append([{'words':[],'factor':[]}])
        for blocks in malformed:
            for fn in (original,compiled):
                with self.assertRaises(ValueError):fn(blocks,10,4)
        for denominator in (0,-1,True,1.5):
            with self.assertRaises(ValueError):compiled(base,denominator,4)


if __name__=='__main__':
    unittest.main()
