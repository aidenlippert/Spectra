"""Exact controls for dual separation, including the off-diagonal convention."""
from fractions import Fraction as F
from itertools import combinations
import unittest
from experiments.marginal_symbolic import mono, add, product, encode, canonical
from research.certificate_scaling.commutator_dual_witness import psd, gram, seed, check


class DualControls(unittest.TestCase):
    def test_psd_zero_pivot_and_negative_directions(self):
        self.assertEqual(psd([[F(1),F(1)],[F(1),F(1)]])['rank'],1)
        self.assertEqual(psd([[F(0),F(0)],[F(0),F(2)]])['rank'],1)
        for bad in ([[0,1],[1,0]], [[1,2],[2,1]], [[-1,0],[0,1]]):
            with self.assertRaises(ValueError):psd([[F(x) for x in row] for row in bad])

    def test_offdiagonal_is_unsymmetrized_moment(self):
        # The physical |+> one-particle density is PSD with off-diagonal 1/2.
        # Using the symmetrized primal coefficient as Y_01 gives 1, a false
        # negative eigenvalue. This catches the factor-of-two audit trap.
        y={((1,i),(0,j)):F(1,2) for i in range(2) for j in range(2)}
        g=gram([mono(((0,0),)),mono(((0,1),))],y)
        self.assertEqual(g,[[F(1,2),F(1,2)],[F(1,2),F(1,2)]])
        self.assertEqual(psd(g)['rank'],1)

    def uniform_witness(self):
        moments=[]
        for k in range(4):
            for inds in combinations(range(4),k):
                w=tuple((1,i) for i in inds)+tuple((0,i) for i in inds)
                moments.append({'word':[list(x) for x in w],'value':str(seed(w,4,2))})
        h=product(mono(((1,0),(0,0))),mono(((1,1),(0,1))))
        return {'modes':4,'particles':2,'hamiltonian':encode(h),'enrich':True,
                'creator_channels':True,'moments':moments}

    def test_uniform_sector_sign_and_full_dual(self):
        self.assertEqual(seed(((1,0),(1,1),(0,0),(0,1)),4,2),F(-1,6))
        rec=check(self.uniform_witness())
        self.assertEqual(F(rec['dual_objective']),F(1,6))

    def test_malformed_and_infeasible_moments_refused(self):
        for mutation in ('box','normalization','duplicate','ideal','enlarged_scope'):
            w=self.uniform_witness()
            if mutation=='box':w['moments'][1]['value']='2'
            elif mutation=='normalization':w['moments'][0]['value']='0'
            elif mutation=='duplicate':w['moments'].append(w['moments'][0].copy())
            elif mutation=='enlarged_scope':w['one_body_steps']=1
            else:w['moments'][1]['value']='1/3'
            with self.assertRaises(ValueError):check(w)


if __name__=='__main__':unittest.main()
