from fractions import Fraction as F
import unittest

import numpy as np

from experiments.marginal_energy_selector import prepare, restricted_problem, exact_certificate, price_candidates
from experiments.marginal_orbit_certificate import average_canonical_polynomial
from experiments.marginal_orbit_quartic import symmetry_group
from experiments.marginal_symbolic import adj, product


class EnergySelectorTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.data=prepare(4)

    def test_pricing_column_matches_exact_averaged_square(self):
        data=self.data;index=max(range(len(data['candidates'])),key=lambda i:len(data['candidates'][i]['words']))
        words=data['candidates'][index]['words'];v=np.array([(i%5-2)/7 for i in range(len(words))])
        p={w:F(i%5-2,7) for i,w in enumerate(words) if i%5!=2}
        exact=average_canonical_polynomial(product(adj(p),p),symmetry_group(4))
        actual=data['candidate_maps'][index]@np.outer(v,v).ravel()
        expected=data['reduced']@np.array([float(exact.get(w,0)) for w in data['rows']])
        self.assertLess(np.max(abs(actual-expected)),1e-12)

    def test_restricted_cubic_problem_exports_valid_interval(self):
        solution=restricted_problem(self.data,[],2)()
        self.assertAlmostEqual(float(self.data['unit']@solution['dual']),1,places=7)
        certificate,receipt=exact_certificate(self.data,[],solution)
        self.assertGreaterEqual(F(receipt['width']),0)
        self.assertLess(F(receipt['width']),F(1,100000))

    def test_reduced_pricing_preserves_full_moment_spectrum(self):
        data=prepare(6,pricing='stabilizer');rng=np.random.default_rng(64)
        for _ in range(3):
            dual=rng.normal(size=len(data['rhs']))
            for i,(candidate,g) in enumerate(zip(data['candidates'],data['candidate_maps'])):
                n=len(candidate['words']);matrix=np.asarray(g.T@dual).reshape(n,n)
                matrix=(matrix+matrix.T)/2;values=[];reconstructed=np.zeros_like(matrix)
                for block in data['pricing_blocks']:
                    if block['candidate']!=i:continue
                    tr=block['transform'];d=tr.shape[1]
                    small=np.asarray(block['map'].T@dual).reshape(d,d);small=(small+small.T)/2
                    reconstructed+=tr@small@tr.T;values.extend(np.linalg.eigvalsh(small))
                self.assertLess(np.max(abs(matrix-reconstructed)),1e-10)
                self.assertTrue(np.allclose(np.linalg.eigvalsh(matrix),sorted(values),atol=1e-10,rtol=0))
            reduced,count=price_candidates(data,dual,1000)
            full,full_count=price_candidates(dict(data,pricing='full'),dual,1000)
            self.assertEqual(count,full_count)
            self.assertEqual({p['candidate'] for p in reduced},{p['candidate'] for p in full})
            for p in reduced:
                i=p['candidate'];v=p['vector']
                self.assertAlmostEqual(float(v@v),1,places=11)
                expected=data['candidate_maps'][i]@np.outer(v,v).ravel()
                self.assertTrue(np.allclose(p['column'],expected,atol=1e-11,rtol=0))
                self.assertAlmostEqual(float(dual@p['column']),p['violation'],places=10)


if __name__=='__main__':unittest.main()
