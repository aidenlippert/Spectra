import unittest
import numpy as np

from experiments.marginal_reynolds import invariant_rows, reynolds_matrix
from experiments.marginal_orbit_quartic import symmetry_group
from experiments.marginal_symbolic import add, canonical, mono, scale, transform, adj


class ReynoldsTests(unittest.TestCase):
    def test_projection_matches_direct_polynomial_average(self):
        modes=4;rows=invariant_rows(modes,4);p,chosen=reynolds_matrix(rows,modes)
        self.assertLess(np.max(abs((p@p-p).toarray())),1e-12)
        for index,w in enumerate(rows):
            polynomial=mono(w)
            if canonical(adj(polynomial))!=polynomial:polynomial=add(polynomial,canonical(adj(polynomial)))
            group=symmetry_group(modes)
            expected=scale(add(*(transform(polynomial,g) for g in group)),1/len(group))
            self.assertLess(np.max(abs(p[:,index].toarray().ravel()-[float(expected.get(r,0)) for r in rows])),1e-12)
        self.assertEqual(len(chosen),np.linalg.matrix_rank(p.toarray()))


if __name__=='__main__':unittest.main()
