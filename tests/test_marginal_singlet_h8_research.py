"""Exact integration/refusal checks for the bounded H8 singlet research CLI."""
import copy,json,unittest
from pathlib import Path
from results.marginal_graded_hubbard8.discovery.singlet_energy_replay import replay,prepare,schur,ldl

ROOT=Path(__file__).resolve().parents[1]/'results/marginal_graded_hubbard8'
class SingletH8ResearchTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.c=json.loads((ROOT/'singlet_energy_first/certificate.json').read_text())
    def test_exact_interval(self):
        r=replay(self.c)
        self.assertTrue(r['accepted']);self.assertEqual(r['lower'],'-6')
        self.assertEqual(r['retained_dimension'],14)
        self.assertLess(r['upper_float'],-2.94)
        self.assertLessEqual(r['unique_action_states'],4096)
    def test_wrong_spin_basis_and_theorem_refused(self):
        c=copy.deepcopy(self.c);c['embedding']['basis'][0][next(iter(c['embedding']['basis'][0]))]=99
        with self.assertRaisesRegex(ValueError,'basis coefficient'):replay(c)
        c=copy.deepcopy(self.c);c['ground_spin_theorem']='SU2-only'
        with self.assertRaisesRegex(ValueError,'theorem premise'):replay(c)
    def test_wrong_hamiltonian_and_response_budget_refused(self):
        c=copy.deepcopy(self.c);c['complement']['U']='5'
        with self.assertRaises(ValueError):replay(c)
        c=copy.deepcopy(self.c);c['response_polynomial']=['1']*9
        with self.assertRaisesRegex(ValueError,'order budget'):replay(c)
    def test_unproved_lower_and_invalid_upper_refused(self):
        c=copy.deepcopy(self.c);c['lower']='-4.3'
        with self.assertRaisesRegex(ValueError,'positive generalized'):replay(c)
        c=copy.deepcopy(self.c);c['independent_upper']={'states':[0],'amplitudes':[1]}
        with self.assertRaises(ValueError):replay(c)
if __name__=='__main__':unittest.main()
