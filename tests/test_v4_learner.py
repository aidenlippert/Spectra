import unittest, numpy as np
from experiments.v4_learner import CARTRegressor, NearestRecordRegressor

class V4Tests(unittest.TestCase):
 def test_unseen_rule_and_threshold(self):
  X=np.arange(32,dtype=float).reshape(-1,1); y=(X[:,0]>13.5).astype(float)
  m=CARTRegressor(max_depth=2,min_leaf=4).fit(X,y); p=m.predict([[2],[20]])
  self.assertEqual(tuple(p),(0.,1.)); self.assertIn("threshold",m.tree)
 def test_roundtrip_and_counters(self):
  m=CARTRegressor(min_leaf=2).fit([[0,0],[1,0],[2,1],[3,1]],[0,0,1,1]); q=CARTRegressor.from_dict(m.to_dict())
  self.assertTrue(np.array_equal(m.predict([[.2,.1],[2.8,.9]]),q.predict([[.2,.1],[2.8,.9]])))
  self.assertGreaterEqual(m.fit_operations,1)
 def test_reject_nan_and_caps(self):
  with self.assertRaises(ValueError): CARTRegressor(max_depth=9)
  with self.assertRaises(ValueError): CARTRegressor().fit([[np.nan]],[1])
  with self.assertRaises(ValueError): CARTRegressor().fit([],[])
  with self.assertRaises(ValueError): CARTRegressor().fit([[1,2]], [1]).predict([[1]])
  with self.assertRaises(ValueError): CARTRegressor.from_dict({"max_depth":1,"min_leaf":1,"n_features":1,"tree":{"value":0,"count":1,"feature":0,"threshold":0,"left":{},"right":{}}})
 def test_nearest_counter(self):
  m=NearestRecordRegressor().fit([[0],[2]],[1,3]); self.assertEqual(tuple(m.predict([[1]])),(1.,)); self.assertEqual(m.predict_operations,2)

if __name__=='__main__': unittest.main()
