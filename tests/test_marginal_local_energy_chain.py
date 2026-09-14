import unittest
from fractions import Fraction as F
from unittest.mock import patch

from experiments.marginal_local_energy_chain import replay, replay_many
from experiments.marginal_local_hubbard_block import replay as replay_block
from experiments.marginal_local_hubbard_block import _terms
from experiments.marginal_transfer_verify import apply_word
from experiments.marginal_polynomial_sos import integer_psd


def singlet_product(sites):
    """Exact product of adjacent valence singlets in the module's CAR order."""
    out = {0: 1}
    for start in range(0, sites, 2):
        local = {1 << (2 * start) | 1 << (2 * (start + 1) + 1): 1,
                 1 << (2 * start + 1) | 1 << (2 * (start + 1)): -1}
        nxt = {}
        for a, ca in out.items():
            for b, cb in local.items():
                nxt[a | b] = nxt.get(a | b, 0) + ca * cb
        out = {k: v for k, v in nxt.items() if v}
    return out


def block(L, U, lower=-10):
    return {'kind': 'local_hubbard_block_v1', 'sites': L, 'U': U, 't': 1,
            'lower': lower, 'upper_vector': singlet_product(L)}


def embedded_terms(N, L, U):
    """Independent symbolic CAR construction of cyclic local operators."""
    v = F(U) * F(L - 1, L)
    terms = []
    for start in range(N):
        for word, coeff in _terms(L, v, 1):
            mapped = tuple((kind, 2 * ((mode // 2 + start) % N) + mode % 2)
                           for kind, mode in word)
            terms.append((mapped, coeff))
        # Centered one-body part and scalar, inserted independently here.
        for local in range(L):
            for spin in (0, 1):
                mode = 2 * ((start + local) % N) + spin
                terms.append((((1, mode), (0, mode)), -v / 2))
        terms.append(((), v * L / 2))
    return terms


def matrix(N, terms):
    result = {}
    for source in range(1 << (2 * N)):
        for word, coeff in terms:
            image = apply_word(word, source)
            if image:
                target, sign = image
                key = (target, source)
                result[key] = result.get(key, F(0)) + coeff * sign
    return {key: value for key, value in result.items() if value}


class LocalEnergyChainTests(unittest.TestCase):
    def test_symbolic_periodic_tiling_identity(self):
        for N, L in ((4, 2), (6, 4)):
            got = matrix(N, embedded_terms(N, L, 4))
            # independently construct centered periodic Hubbard operator
            expected_terms = []
            for i in range(N):
                expected_terms.append((((1, 2*i), (0, 2*i), (1, 2*i+1), (0, 2*i+1)), F(4)))
                expected_terms.append((((1, 2*i), (0, 2*i)), F(-2)))
                expected_terms.append((((1, 2*i+1), (0, 2*i+1)), F(-2)))
                j = (i + 1) % N
                for s in (0, 1):
                    a, b = 2*i+s, 2*j+s
                    expected_terms += [(((1,a),(0,b)), -1), (((1,b),(0,a)), -1)]
            # The centered interaction contains the scalar +U/2 per site.
            expected_terms.append(((), F(2 * N)))
            expected = matrix(N, expected_terms)
            scale = L - 1
            self.assertEqual({key: value / scale for key, value in got.items()}, expected)

    def test_wrap_bond_two_t_norm_bound_is_psd(self):
        N = 2
        terms = []
        for spin in (0,1):
            terms += [(((1,spin),(0,2+spin)),-1),(((1,2+spin),(0,spin)),-1)]
        T = matrix(N, terms)
        for sign in (1, -1):
            dim = 1 << (2 * N)
            integer = [[2 * (i == j) + sign * T.get((i, j), 0)
                        for j in range(dim)] for i in range(dim)]
            integer = [[int(x) for x in row] for row in integer]
            self.assertTrue(integer_psd(integer)['positive_semidefinite'])

    def test_global_tiling_remainder_and_overlap(self):
        cert = {'kind': 'local_energy_chain_v1', 'sites': 10, 'block_length': 4,
                'U': 4, 't': 1, 'blocks': [block(4, 4), block(2, 4)],
                'overlap': block(4, F(3), lower=-10)}
        got = replay(cert)
        self.assertEqual(got['full_blocks'], 2)
        self.assertEqual(got['remainder_sites'], 2)
        self.assertEqual(got['overlap_lower'], '-106/3')
        self.assertEqual(got['upper'], '0')
        huge = dict(cert, sites=10**9)
        huge['blocks'] = [block(4, 4)]
        huge.pop('overlap')
        self.assertEqual(replay(huge)['full_blocks'], 250000000)

    def test_structural_rejections(self):
        base = {'kind': 'local_energy_chain_v1', 'sites': 10, 'block_length': 4,
                'U': 4, 't': 1, 'blocks': [block(4, 4), block(2, 4)]}
        for bad in ({**base, 'blocks': [block(4, 4), block(4, 4)]},
                    {**base, 'blocks': [block(4, 4)]},
                    {**base, 'sites': 9},
                    {**base, 'overlap': block(4, 4)},
                    {**base, 'blocks': [{**block(4,4),'upper_vector':{'0':1}},block(2,4)]},
                    {**base, 'blocks': [block(4, 4), block(2, 4, lower=1)]}):
            with self.assertRaises(ValueError):
                replay(bad)

    def test_batch_cache_recomputes_unique_local_proofs_and_rejects_bad_bounds(self):
        c = {'kind':'local_energy_chain_v1','sites':8,'block_length':2,
             'U':4,'t':1,'blocks':[block(2,4)]}
        with patch('experiments.marginal_local_energy_chain.replay_block',wraps=replay_block) as checked:
            result = replay_many([c,{**c,'sites':1000000}])
            self.assertEqual(checked.call_count,1)
            self.assertEqual(result['unique_local_certificates'],1)
        with self.assertRaises(ValueError):
            replay_many([c,{**c,'blocks':[block(2,4,lower=0)]}])
        with self.assertRaises(ValueError):replay_many([])

    def test_weighted_window_and_distinct_tiling_size(self):
        overlap = {**block(4,3,lower='-2.04052'),
                   'onsite_profile':['1/2','11/2','11/2','1/2'],
                   'hopping_profile':['3/4','3/2','3/4']}
        c = {'kind':'local_energy_chain_v1','sites':1000000,'block_length':2,
             'U':4,'t':1,'blocks':[block(2,4)],'overlap':overlap}
        result = replay(c)
        self.assertEqual(F(result['overlap_lower']),F(1000000,3)*F('-2.04052')-2)
        bad_tile = {**block(4,4),'onsite_profile':[1,7,7,1]}
        with self.assertRaisesRegex(ValueError,'uniform chain'):
            replay({**c,'block_length':4,'blocks':[bad_tile]})
        # Independent full CAR equality for nonuniform cyclic windows.
        N,L = 6,4
        terms = []
        for start in range(N):
            for word,a in _terms(L,3,1,list(map(F,overlap['onsite_profile'])),list(map(F,overlap['hopping_profile']))):
                terms.append((tuple((kind,2*((mode//2+start)%N)+mode%2) for kind,mode in word),a))
            for i,U in enumerate(map(F,overlap['onsite_profile'])):
                terms.append(((),U/2))
                for spin in (0,1):
                    mode=2*((start+i)%N)+spin
                    terms.append((((1,mode),(0,mode)),-U/2))
        self.assertEqual(matrix(N,terms),matrix(N,embedded_terms(N,L,4)))


if __name__ == '__main__':
    unittest.main()
