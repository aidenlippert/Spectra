import unittest
from fractions import Fraction as F
from unittest.mock import patch
from experiments.marginal_local_hubbard_block import replay, sector_matrices, _reflection, _actions
from experiments.marginal_hubbard_polynomial import build
from experiments.marginal_symbolic import decode
from experiments.marginal_transfer_verify import apply_word


def certificate(**updates):
    result = {'kind':'local_hubbard_block_v1','sites':2,'U':'4','t':'1','lower':'-0.82844',
              'upper_vector':{'9':1,'6':-1}}
    result.update(updates)
    return result


class LocalBlockTests(unittest.TestCase):
    def test_complete_basis_and_original_car_congruence(self):
        for L,U,t in [(2,F(4),F(1)),(4,F(3,2),F(2,3))]:
            data = sector_matrices(L,U,t)
            self.assertEqual(sum(len(x[1]) for x in data.values()),4**L)
            h = decode(build(L,U,t)['hamiltonian'],2*L,4)
            physical = {}
            for state in range(4**L):
                image = {state:-U*F(state.bit_count()-L,2)}
                for word,a in h.items():
                    out = apply_word(word,state)
                    if out:image[out[0]] = image.get(out[0],0)+a*out[1]
                physical[state] = {s:a for s,a in image.items() if a}
            self.assertEqual(_actions(L,U,t),physical)
            all_columns = []
            for (up,down,parity),(states,matrix,columns) in data.items():
                all_columns.extend(columns)
                for i,v in enumerate(columns):
                    reflected = {}
                    for s,a in v.items():
                        target,sign = _reflection(s,L)
                        reflected[target] = a*sign
                    self.assertEqual(reflected,{s:parity*a for s,a in v.items()})
                    for j,w in enumerate(columns):
                        actual = sum(a*b*physical[s].get(t,0) for s,a in w.items() for t,b in v.items())
                        self.assertEqual(matrix[i][j],actual)
            for i,v in enumerate(all_columns):
                for j,w in enumerate(all_columns):
                    gram = sum(a*w.get(s,0) for s,a in v.items())
                    self.assertEqual(gram,sum(a*a for a in v.values()) if i==j else 0)

    def test_tight_threshold_uses_amplitude_gram_not_state_labels(self):
        result = replay(certificate())
        self.assertEqual(result['sum_dimensions'],16)
        self.assertEqual(result['upper_quotient'],'0')
        self.assertEqual(result['upper_norm'],'2')
        # Exact two-site ground energy is2-2sqrt(2), between these endpoints.
        with self.assertRaises(ValueError):replay(certificate(lower='-0.82842'))
        with self.assertRaises(ValueError):replay(certificate(lower='0'))

    def test_six_site_local_coverage_respects_existing_psd_cap(self):
        data = sector_matrices(6,4,1)
        self.assertEqual(sum(len(item[1]) for item in data.values()),4096)
        self.assertEqual(len(data),94)
        self.assertEqual(max(len(item[1]) for item in data.values()),200)

    def test_weighted_profiles_and_centered_upper(self):
        onsite = [F(1,2),F(11,2),F(11,2),F(1,2)]
        hopping = [F(3,4),F(3,2),F(3,4)]
        actions = _actions(4,3,1,onsite,hopping)
        # Independently construct physical hopping and centered onsite values.
        for state in range(256):
            diagonal = sum(onsite[i]*(((state>>(2*i))&3)==3)
                           -F(onsite[i],2)*(((state>>(2*i))&3).bit_count()-1)
                           for i in range(4))
            expected = {state:diagonal}
            for i,t in enumerate(hopping):
                for spin in (0,1):
                    x,y = 2*i+spin,2*(i+1)+spin
                    for word in (((1,x),(0,y)),((1,y),(0,x))):
                        out = apply_word(word,state)
                        if out:expected[out[0]] = expected.get(out[0],0)-t*out[1]
            self.assertEqual(actions[state],{s:a for s,a in expected.items() if a})
        c = certificate(sites=4,U=3,onsite_profile=onsite,hopping_profile=hopping,
                        lower='-2.04052',upper_vector={'195':1})
        result = replay(c)
        # Both edge sites doubled, both middle sites empty: physical U-energy1,
        # but the centered profile operator has diagonal energy6.
        self.assertEqual(result['upper_quotient'],'6')
        with self.assertRaises(ValueError):replay({**c,'lower':'-2.04051'})
        for update in [dict(onsite_profile=[1,5,5,1.0]),
                       dict(onsite_profile=[1,5,5,2]),
                       dict(onsite_profile=[1,5,4,2]),
                       dict(hopping_profile=[-1,5,-1]),
                       dict(hopping_profile=[1,1])]:
            with self.assertRaises(ValueError):replay({**c,**update})

    def test_noncommuting_reflection_is_rejected(self):
        # Swap just states0/1, which does not commute with the centered model.
        def broken(state,L):return (1-state,1) if state in (0,1) else (state,1)
        with patch('experiments.marginal_local_hubbard_block._reflection',broken):
            with self.assertRaisesRegex(ValueError,'commute'):replay(certificate())

    def test_refusals(self):
        for updates in [dict(kind='bad'),dict(sites=8),dict(sites=True),dict(U=4.0),
                        dict(t=-1),dict(lower=-1.0),dict(lower='1/10000000000'),
                        dict(upper_vector={}),dict(upper_vector={'9':0}),
                        dict(upper_vector={'0':1}),dict(upper_vector={'16':1}),
                        dict(upper_vector={'9':1.0}),dict(upper_vector={'9':10**13}),
                        dict(upper_vector={'9':1,9:1})]:
            with self.assertRaises(ValueError):replay(certificate(**updates))

if __name__=='__main__':unittest.main()
