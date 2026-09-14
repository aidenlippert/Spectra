import unittest
from fractions import Fraction as F
from experiments.v7_prediction_requests import request_prediction, _feature
from experiments.v7_identifiability import verify_compatible_models


class PredictionRequestTests(unittest.TestCase):
    def setUp(self):
        self.X = [[F(1), F(0), F(1)], [F(1), F(0), F(1)]]
        self.y = [F(-49), F(-49)]
        self.theta = [F(1), F(0), F(-50)]

    def test_ambiguous_first_query_blocks_horizon(self):
        out = request_prediction({'X': self.X, 'y': self.y, 'theta0': self.theta, 'order': 1, 'inputs': 1},
                                 [F(0)], [], [F(1)], F(1), F(0), order=1, input_width=1)
        self.assertEqual(out['status'], 'no_guarantee')
        w = out['witness']
        self.assertTrue(verify_compatible_models(self.X, self.y, w['query'],
                                                  w['theta_plus'], w['theta_minus'], F(0), F(2)))

    def test_rowspace_query_is_supported(self):
        X = [[F(1), F(0)], [F(0), F(1)]]
        out = request_prediction({'X': X, 'y': [F(2), F(3)], 'theta0': [F(2), F(3)], 'order': 1, 'inputs': 0},
                                 [F(0)], [], [], F(1), F(0), order=1, input_width=0)
        self.assertEqual(out['status'], 'guaranteed')

    def test_actual_v6_feature_order_and_history(self):
        self.assertEqual(_feature([2,3],[[7,11]],[13,17],2,2),
                         list(map(F,[3,2,13,17,7,11,1])))
        self.assertNotEqual(_feature([2,3],[[7,11]],[13,17],2,2),
                            _feature([2,4],[[7,11]],[13,17],2,2))

    def test_rank_deficient_same_input_and_wide_support(self):
        data=dict(X=[[1,50,1],[1,50,1]],y=[1,1],theta0=[1,0,0],order=1,inputs=1)
        a=request_prediction(data,[1],[],[50],F(1,10),F(1,20),witness=[0,1,-50])
        self.assertEqual(a['status'],'guaranteed')
        self.assertEqual(F(a['support_error_bound']),F(1,20))
        a=request_prediction(data,[1],[],[50],F(1,100),F(1,20),witness=[0,1,-50])
        self.assertEqual(a['status'],'no_guarantee')
        a=request_prediction(data,[1],[],[75],F(1,10),F(1,20),witness=[0,1,-50],horizon=16)
        self.assertEqual(a['status'],'no_guarantee')
        self.assertEqual(a['horizon'],16)
        self.assertIsNotNone(a['witness'])

    def test_actual_measured_record_witness(self):
        from experiments.v7_measured_gate import run
        r=run()
        self.assertTrue(r['exact_witness_verified'])
        self.assertFalse(r['future_response_outputs_read'])
        self.assertEqual(r['training_rows'],897)
        self.assertGreater(F(r['prediction_disagreement']),2*F(r['tolerance']))

    def test_no_future_outputs_in_api_and_short_history_abstains(self):
        out = request_prediction({'X': [[F(1), F(0), F(1)]], 'y': [F(0)],
                                  'theta0': [F(0), F(0), F(0)], 'order': 1, 'inputs': 2},
                                 [], [], [F(1), F(0)], F(1), F(0), order=1)
        self.assertEqual(out['status'], 'no_guarantee')


if __name__ == '__main__':
    unittest.main()
