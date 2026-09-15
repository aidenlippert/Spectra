import unittest
import copy,json
from pathlib import Path
from fractions import Fraction as F
from experiments.marginal_window_family_bound import affine_rayleigh, local_certificate, replay
from experiments.marginal_local_hubbard_block import _actions


class WindowFamilyBoundTests(unittest.TestCase):
    def test_affine_rayleigh_matches_actual_profile_actions(self):
        vector={23:7,29:-2,39:3,45:5}
        norm,(c,alpha,beta)=affine_rayleigh(vector)
        for a,b in [(F(0),F(0)),(F(6),F(3,2)),(F(1,2),F(3,4)),(F(7,5),F(2,3))]:
            actions=_actions(4,3,1,[a,6-a,6-a,a],[b,3-2*b,b])
            energy=sum(value*coefficient*vector.get(target,0) for source,value in vector.items()
                       for target,coefficient in actions[source].items())/F(norm)
            self.assertEqual(c+a*alpha+b*beta,energy)

    def test_exact_rectangle_upper(self):
        # One doubled edge, opposite empty edge, two singly occupied middle
        # sites gives centered expectation exactly a for every b.
        norm,coefficients=affine_rayleigh({'23':1})
        self.assertEqual(norm,1)
        self.assertEqual(coefficients,(0,1,0))
        certificate={'kind':'hubbard_four_window_family_bound_v1','a':'1/2','b':'3/4',
                     'lower':'-2.04052','upper_vector':{'23':1},'target_family_upper':'6'}
        result=replay(certificate)
        self.assertEqual(result['maximum_local_minimum_upper'],'6')
        with self.assertRaisesRegex(ValueError,'family upper'):
            replay({**certificate,'target_family_upper':'5'})
        with self.assertRaises(ValueError):replay({**certificate,'lower':'-2.04051'})

    def test_refusal_gates(self):
        for vector in [{},{'23':0},{'0':1},{'23':1.0},{'23':10**13},{'23':1,23:2}]:
            with self.assertRaises(ValueError):affine_rayleigh(vector)
        for a,b in [(-1,1),(7,1),(1,2),(True,1),(1,0.5)]:
            with self.assertRaises(ValueError):local_certificate(a,b,-10)
        with self.assertRaises(ValueError):replay({'kind':'wrong'})

    def test_sharp_fixture_and_actual_witness_tamper(self):
        path=Path(__file__).resolve().parents[1]/'results/marginal_graded_hubbard8/weighted_window_family/certificate.json'
        certificate=json.loads(path.read_text())
        result=replay(certificate)
        self.assertLess(F(result['family_bracket_width']),F(1,10**9))
        self.assertEqual(result['affine_expectation'][1:],['0','0'])
        consistency=result['overlap_consistency']
        self.assertTrue(consistency['three_site_reductions_equal'])
        self.assertLess(F(consistency['three_site_purity']),1)
        self.assertTrue(consistency['no_five_site_translation_invariant_extension'])
        self.assertEqual(result['arbitrary_three_site_boundary_correction_upper'],result['maximum_local_minimum_upper'])
        bad=copy.deepcopy(certificate)
        key=max(bad['upper_vector'],key=lambda s:abs(bad['upper_vector'][s]))
        bad['upper_vector'][key]//=2
        with self.assertRaisesRegex(ValueError,'family upper'):replay(bad)

if __name__=='__main__':unittest.main()
