"""Independent CAR action and reference-witness checks for integer streaming."""
from fractions import Fraction as F
from itertools import combinations
import unittest
from experiments.marginal_symbolic import encode,number_shift
from experiments.marginal_transfer_verify import apply_word
from experiments.marginal_determinant_tree import DeterminantOracle
from research.certificate_scaling.streaming_reference_upper import compile_term,upper

class StreamingReferenceControls(unittest.TestCase):
    def test_compiled_actions_match_ladder_application(self):
        for k in range(3):
            for left in combinations(range(4),k):
                for right in combinations(range(4),k):
                    word=tuple((1,i) for i in left)+tuple((0,i) for i in right)
                    req,occ,flip,parity,c=compile_term(word,1)
                    for state in range(16):
                        actual=apply_word(word,state)
                        predicted=(state^flip,c*(-1 if (state&parity).bit_count()%2 else 1)) if state&req==occ else None
                        self.assertEqual(predicted,actual)

    def test_large_diagonal_witness_without_action_cache(self):
        cert={'modes':16,'particles':8,'hamiltonian':encode(number_shift(16,0))}
        states=[sum(1<<i for i in occupied) for occupied in combinations(range(16),8)]
        witness={'states':states,'amplitudes':[1]*len(states)}
        with self.assertRaises(ValueError):DeterminantOracle(cert).upper(witness)
        value,rec=upper(cert,witness)
        self.assertEqual(value,F(8));self.assertEqual(rec['action_cache_states'],0)
        self.assertEqual(rec['witness_states'],12870)

    def test_validation_refusals(self):
        cert={'modes':4,'particles':2,'hamiltonian':encode(number_shift(4,0))}
        for witness in ({'states':[3,3],'amplitudes':[1,1]}, {'states':[3],'amplitudes':[0]},
                        {'states':[1],'amplitudes':[1]}, {'states':[3],'amplitudes':[1.]},
                        {'states':[3]*65537,'amplitudes':[1]*65537}):
            with self.assertRaises(ValueError):upper(cert,witness)

if __name__=='__main__':unittest.main()
