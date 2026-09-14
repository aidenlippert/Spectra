from fractions import Fraction as F
import unittest

import numpy as np

from experiments.marginal_adaptive import assemble_and_solve
from experiments.marginal_coefficient import dictionaries, export
from experiments.marginal_symbolic import canonical, number_shift, validate_word, word_product


class AdaptiveTests(unittest.TestCase):
    def test_cached_boolean_labels_do_not_leak_into_validated_words(self):
        word_product(((True,0),(False,0)),())
        p=canonical({((1,0),(0,0)):1})
        for w in p:validate_word(w,4,4)

    def test_residual_penalty_keeps_certified_number_bound(self):
        h=number_shift(4,0);blocks=dictionaries(4,'quadratic')
        proposal,data=assemble_and_solve(h,4,2,blocks,residual_penalty=True)
        _,receipt=export(h,4,2,blocks,proposal)
        self.assertGreater(F(receipt['lower']),F(1999,1000))
        self.assertLessEqual(F(receipt['lower']),2)
        self.assertAlmostEqual(proposal['numerical_objective'],2,places=6)

    def test_rank_one_polynomial_block_exports_exactly(self):
        h=number_shift(4,0);blocks=dictionaries(4,'quadratic')
        blocks.append({'name':'polynomial','words':[(),((1,0),(0,0))],
                       'basis_transform':[[1],[1]]})
        proposal,data=assemble_and_solve(h,4,2,blocks)
        _,receipt=export(h,4,2,blocks,proposal)
        self.assertGreater(F(receipt['lower']),F(1999,1000))
        self.assertEqual(proposal['gram_dimensions'][-1],1)

    def test_primal_and_moment_proposals_obey_number_problem(self):
        h=number_shift(4,0);blocks=dictionaries(4,'quadratic')
        proposal,data=assemble_and_solve(h,4,2,blocks)
        cert,receipt=export(h,4,2,blocks,proposal)
        self.assertGreater(F(receipt['lower']),F(1999,1000))
        y=data['dual']
        self.assertAlmostEqual(y[data['rows'].index(())],1,places=7)
        self.assertAlmostEqual(float(data['rhs']@y),2,places=7)
        self.assertLess(np.max(np.abs(data['multiplier_map'].T@y)),1e-7)
        for g,block in zip(data['maps'],blocks):
            matrix=(g.T@y).reshape(len(block['words']),len(block['words']))
            # The one-representative Hermitian coefficient map need not return
            # a symmetric raw matrix; only its symmetric part enters a Gram.
            self.assertGreaterEqual(np.linalg.eigvalsh((matrix+matrix.T)/2)[0],-1e-7)

    def test_rejects_invalid_symmetry_or_degree(self):
        h=number_shift(4,0);blocks=dictionaries(4,'quadratic')
        with self.assertRaises(ValueError):assemble_and_solve(h,4,2,blocks,degree=5)
        with self.assertRaises(ValueError):assemble_and_solve(h,4,2,blocks,groups=[[0,1],[1,2,3]])
        with self.assertRaises(ValueError):assemble_and_solve(h,4,2,blocks,groups=[[0],[1],[2],[3]])
        blocks.append({'words':[()], 'basis_transform':[[float('nan')]]})
        with self.assertRaises(ValueError):assemble_and_solve(h,4,2,blocks)

    def test_exact_binary_symmetry_filter(self):
        h=number_shift(4,0);blocks=[{'name':'constant','words':[()]}]
        proposal,data=assemble_and_solve(h,4,2,blocks,parity_masks=[1])
        _,receipt=export(h,4,2,blocks,proposal)
        self.assertGreater(F(receipt['lower']),F(1999,1000))
        self.assertTrue(all(sum(i==0 for c,i in w)%2==0 for w in data['rows']))
        hopping={((1,0),(0,1)):F(1),((1,1),(0,0)):F(1)}
        with self.assertRaises(ValueError):assemble_and_solve(hopping,4,2,blocks,parity_masks=[1])
        with self.assertRaises(ValueError):assemble_and_solve(h,4,2,blocks,parity_masks=[True])
        mixed=[{'name':'mixed','words':[((0,0),),((0,1),)]}]
        with self.assertRaises(ValueError):assemble_and_solve(h,4,2,mixed,parity_masks=[1])


if __name__=='__main__':unittest.main()
