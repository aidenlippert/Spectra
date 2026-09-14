"""Compare every proposed map coefficient with independent exact CAR products."""
import unittest
from fractions import Fraction as F
import numpy as np
from experiments.marginal_symbolic import mono,add,scale
from research.certificate_scaling.commutator_dictionary import gram_columns
from research.certificate_scaling.polynomial_gram_contraction import prepare,contract,whitening_transform


class ContractionControls(unittest.TestCase):
    def test_transformed_map_matches_independent_exact_polynomials(self):
        p=add(mono(((0,0),),F(1,3)),mono(((1,1),(0,0),(0,2)),F(2,7)))
        q=add(mono(((0,2),),F(-3,5)),mono(((1,0),(0,1),(0,2)),F(1,11)))
        gs=[{'polynomials':[p,q]}];blocks,words,_=prepare(gs)
        rows=sorted(words,key=lambda w:(len(w),w));lookup={w:i for i,w in enumerate(rows)}
        W=np.array([[2.,1.],[0.,1.]])
        actual,indices,_=contract(blocks[0],[p,q],lookup,transform=W)
        exact,expected_indices,_=gram_columns([add(scale(p,2),q),q])
        expected=np.array([[float(c.get(w,0)) for c in exact] for w in rows])
        np.testing.assert_allclose(actual.toarray(),expected,atol=1e-15,rtol=1e-14)
        self.assertEqual(indices,expected_indices)
        with self.assertRaises(ValueError):contract(blocks[0],[p,q],lookup,transform=np.ones((1,2)))

    def test_whitening_preserves_dimension_and_refuses_rank_truncation(self):
        p=mono(((0,0),));q=add(p,mono(((0,1),),F(1,1000)))
        W,r=whitening_transform([p,q])
        C=np.array([[1.,1.],[0.,.001]])
        np.testing.assert_allclose((C@W.T).T@(C@W.T),np.eye(2),atol=1e-10)
        self.assertEqual(W.shape,(2,2));self.assertGreater(r['estimated_condition_before'],1000)
        with self.assertRaises(ValueError):whitening_transform([p,p])

    def test_whitening_keeps_disjoint_support_zeros_exact(self):
        p=mono(((0,0),));q=add(p,mono(((0,1),),F(1,1000)))
        r=mono(((0,2),));t=add(r,mono(((0,3),),F(1,1000)))
        W,stats=whitening_transform([p,q,r,t])
        self.assertEqual(stats['support_component_dimensions'],[2,2])
        self.assertTrue(np.all(W[:2,2:]==0));self.assertTrue(np.all(W[2:,:2]==0))
        blocks,words,_=prepare([{'polynomials':[p,q,r,t]}])
        lookup={w:i for i,w in enumerate(sorted(words))}
        matrix,indices,_=contract(blocks[0],[p,q,r,t],lookup,transform=W)
        # Coefficient-orthogonal components are not physically independent:
        # cross-component Gram variables MUST remain available.
        self.assertGreater(matrix[:,indices.index(2)].nnz,0) # Q[0,2]

    def test_full_map_diagonal_offdiagonal_and_degree_six(self):
        p=add(mono(((0,0),),F(1,3)),mono(((1,1),(0,0),(0,2)),F(2,7)))
        q=add(mono(((0,2),),F(-3,5)),mono(((1,0),(0,1),(0,2)),F(1,11)))
        groups=[{'polynomials':[p,q]}];blocks,words,stats=prepare(groups)
        rows=sorted(words,key=lambda w:(len(w),w));lookup={w:i for i,w in enumerate(rows)}
        actual,indices,_=contract(blocks[0],[p,q],lookup,batch=1)
        exact,expected_indices,_=gram_columns([p,q])
        expected=np.array([[float(c.get(w,0)) for c in exact] for w in rows])
        np.testing.assert_allclose(actual.toarray(),expected,atol=1e-15,rtol=1e-14)
        self.assertEqual(indices,expected_indices)
        self.assertTrue(any(len(w)==6 for w in words))
        self.assertEqual(stats['monomial_word_pairs'],16)
        with self.assertRaises(ValueError):contract(blocks[0],[p,q],lookup,batch=0)


if __name__=='__main__':unittest.main()
