import copy
from fractions import Fraction as F
import unittest

from research.side_routes_20260913.combine_bounds import combine
from research.side_routes_20260913.finite_range import bounds


class CombinedCertificateTests(unittest.TestCase):
    def payloads(self):
        model = {'modes':3,'particles':1,'hopping':['1','1'],'interaction':['0','0'],'fields':['0','1','3']}
        result = []
        for sites in (['1','1','1'],['2','1','1']):
            amp = {'sites':sites,'pairs':[]}
            result.append({'model':model,'amplitude':amp,'claim':bounds(model,amp)})
        return result

    def test_two_distinct_states_give_strictly_narrower_combination(self):
        payloads = self.payloads(); result = combine(payloads)
        self.assertEqual(F(result['lower']),-1)
        self.assertEqual(F(result['upper']),F(-1,3))
        self.assertEqual(F(result['width']),F(2,3))
        self.assertTrue(all(F(result['width'])<F(p['claim']['width']) for p in payloads))

    def test_rejects_model_mismatch_and_tampered_winner(self):
        payloads = self.payloads(); payloads[1]['model'] = copy.deepcopy(payloads[1]['model'])
        payloads[1]['model']['fields'][0]='100'
        with self.assertRaises(ValueError): combine(payloads)
        payloads = self.payloads(); payloads[0]['claim']['lower']='100'
        with self.assertRaises(ValueError): combine(payloads)


class LowerSearchTests(unittest.TestCase):
    def test_lower_objective_preserves_or_improves_initial_bound(self):
        from research.side_routes_20260913.range_campaign import search
        result = search(('repulsive',1,'lower',17))
        initial = next(p for p in result['candidates'] if p['candidate']=='initial')
        self.assertGreaterEqual(F(result['winner']['claim']['lower']),F(initial['claim']['lower']))
        self.assertGreater(result['numerical_objective_calls'],10)
        with self.assertRaises(ValueError): search(('repulsive',1,'unsupported',17))


if __name__ == '__main__': unittest.main()
