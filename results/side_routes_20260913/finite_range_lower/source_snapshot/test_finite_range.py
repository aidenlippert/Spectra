import copy
from fractions import Fraction as F
from itertools import combinations
import unittest

from research.side_routes_20260913.finite_range import bounds, compile_factors, replay
from research.side_routes_20260913.positive_chain import bounds as nearest_bounds
from research.side_routes_20260913.test_positive_chain import independent


class FiniteRangeTests(unittest.TestCase):
    def test_every_sector_against_independent_CAR_and_full_amplitudes(self):
        m = 6
        model = {'modes': m, 'particles': 0, 'hopping': ['1','2/3','3/5','5/7','0'],
                 'interaction': ['1/3','-2/5','2/7','3/4','-1/2'],
                 'fields': ['1/7','-1/3','2/5','-3/7','4/9','-5/11'], 'offset': '2/13'}
        sites = ['1','3/2','4/3','5/4','6/5','7/6']
        pairs = [['2/3','3/4','4/5','5/6','6/7'], ['3/2','4/3','5/4','6/5'], ['2','3/2','4/3']]
        for radius in range(4):
            amp = {'sites': sites, 'pairs': pairs[:radius]}
            for n in range(m + 1):
                model['particles'] = n
                states, H, *_ = independent(model, {'sites': sites, 'bonds': ['1'] * (m - 1)})
                psi = []
                for state in states:
                    p = F(1)
                    for i in range(m): p *= F(sites[i]) ** ((state >> i) & 1)
                    for d, row in enumerate(amp['pairs'], 1):
                        for i, x in enumerate(row): p *= F(x) ** (((state >> i) & 1) * ((state >> (i+d)) & 1))
                    psi.append(p)
                local = [sum(H[i][j]*psi[j] for j in range(len(states))) / psi[i] for i in range(len(states))]
                Z = sum(p*p for p in psi)
                upper = sum(p*p*e for p,e in zip(psi,local)) / Z
                result = bounds(model, amp)
                self.assertEqual(F(result['lower']), min(local))
                self.assertEqual(F(result['upper']), upper)
                self.assertEqual(F(result['partition']), Z)
                ending, _ = compile_factors(model, amp)
                for state, expected in zip(states, local):
                    tail = 0; actual = F(model['offset'])
                    for i in range(m):
                        tail = (tail << 1) | ((state >> i) & 1)
                        actual += sum(table[tail & ((1 << width)-1)] for width,table in ending[i])
                    self.assertEqual(actual, expected)

    def test_range_one_matches_original_and_padding_preserves_bound(self):
        m = 8
        model = {'modes': m, 'particles': 4, 'hopping': ['1']*7,
                 'interaction': ['1/2']*7, 'fields': ['0']*8}
        sites = ['1','3/2']*4; pairs = [['3/4']*7]
        original = nearest_bounds(model, {'sites':sites, 'bonds':pairs[0]})
        for radius in (1,2,3):
            amp = {'sites':sites, 'pairs':pairs+[['1']*(m-d) for d in range(2,radius+1)]}
            result = bounds(model, amp)
            for key in ('lower','upper','width','partition'): self.assertEqual(original[key],result[key])

    def test_numerical_proposals_match_exact_dyadic_witnesses(self):
        import numpy as np
        from research.side_routes_20260913.range_discovery import NumericalChain,rational_amplitude
        from research.side_routes_20260913.chain_campaign import case
        for name in ('repulsive','site_disorder','bond_disorder'):
            model = case(name, 8)
            for radius in range(4):
                parameters = np.array([-.35,.15,-.1,.06,.12,-.04][:radius+3])
                amp, weights = rational_amplitude(model,radius,parameters)
                actual_parameters = np.log([float(F(w)) for w in weights])
                lo, hi = NumericalChain(model,radius).interval(actual_parameters)
                result = bounds(model,amp)
                self.assertAlmostEqual(lo,result['lower_float'],places=10)
                self.assertAlmostEqual(hi,result['upper_float'],places=10)

    def test_refuses_invalid_witnesses_and_false_claims(self):
        model = {'modes':4,'particles':2,'hopping':['1']*3,'interaction':['1']*3,'fields':['0']*4}
        amp = {'sites':['1']*4,'pairs':[['1']*3]}
        payload = {'model':model,'amplitude':amp,'claim':bounds(model,amp)}
        replay(payload)
        bad = copy.deepcopy(payload); bad['claim']['lower']='1000'
        with self.assertRaises(ValueError): replay(bad)
        for pairs in ([['1']*2], [['1','0','1']], [[1.,1,1]], [['1']*3,['1']*2,['1'],[]]):
            with self.assertRaises(ValueError): bounds(model,{'sites':amp['sites'],'pairs':pairs})
        bad = copy.deepcopy(model); bad['hopping'][0]='-1'
        with self.assertRaises(ValueError): bounds(bad,amp)


if __name__ == '__main__': unittest.main()
