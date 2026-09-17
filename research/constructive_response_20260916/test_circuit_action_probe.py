import unittest
import numpy as np
from research.constructive_response_20260916.circuit_action_probe import setup


class CircuitActionTests(unittest.TestCase):
    def test_matches_expanded_action_for_vector_and_columns(self):
        labels,h,gates,action=setup(3)
        expanded=h
        for g in gates: expanded=(g.T@expanded@g).tocsr()
        v=np.sin(np.arange(len(labels))+0.17)
        block=np.column_stack([v,np.cos(np.arange(len(labels)))])
        np.testing.assert_allclose(action(v),expanded@v,rtol=3e-13,atol=3e-13)
        np.testing.assert_allclose(action(block),expanded@block,rtol=3e-13,atol=3e-13)


if __name__=='__main__':
    unittest.main()
