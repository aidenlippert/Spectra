import unittest
from experiments.marginal_valence_states import valence_states,dimer_singlet_witness,spin_raise_residual
from experiments.marginal_spin_reduction import SpinZeroOracle
from experiments.marginal_spin_constructor import spin_states

class ValenceStatesTests(unittest.TestCase):
    def test_direct_states_equal_single_occupation_projector(self):
        for modes,count in ((4,2),(8,6),(12,20)):
            p=valence_states(modes,modes//2);self.assertEqual(len(p),count)
            o=SpinZeroOracle({'modes':modes,'particles':modes//2,'hamiltonian':[]})
            expected=[s for s in spin_states(o) if all(((s>>(2*i))&3) in (1,2) for i in range(modes//2))]
            self.assertEqual(p,expected)
    def test_dimer_witness_is_exact_nonzero_singlet(self):
        for modes in (4,8,12):
            p=valence_states(modes,modes//2);w=dimer_singlet_witness(modes,modes//2)
            self.assertEqual(len(w['states']),2**(modes//4));self.assertTrue(set(w['states'])<=set(p))
            self.assertEqual(spin_raise_residual(w,modes),{})
            self.assertEqual(sum(a*a for a in w['amplitudes']),len(w['states']))
            for s in w['states']: self.assertEqual(s.bit_count(),modes//2)
    def test_limits_fail_explicitly_instead_of_truncating(self):
        for modes,particles in ((True,2),(4,True),(0,0),(6,3),(12,4),(16,8),(68,34)):
            with self.assertRaises(ValueError): valence_states(modes,particles)
            with self.assertRaises(ValueError): dimer_singlet_witness(modes,particles)
        for cap in (True,0,19,33):
            with self.assertRaises(ValueError): valence_states(12,6,cap)

if __name__=='__main__':unittest.main()
