import unittest
import json
from fractions import Fraction
from experiments.v2_policy import certificate,verify
class PolicyTests(unittest.TestCase):
 def test_h1_reference_values(self):
  self.assertEqual(certificate("none",1)["value"],"7/16"); self.assertEqual(certificate("A",1)["value"],"3/8"); self.assertEqual(certificate("B",1)["value"],"3/8"); self.assertEqual(certificate("AB",1)["value"],"1/4")
 def test_all_actions_and_tamper(self):
  c=certificate("none",1); self.assertTrue(verify(c)); c["value"]="0"; self.assertFalse(verify(c))
 def test_horizon(self):
  c=certificate("none",3); self.assertTrue(verify(c)); self.assertLessEqual(Fraction(c["value"]),Fraction(7,16)); c["nodes"][0]["value"]="0"; self.assertFalse(verify(c))
 def test_json_and_structural_corruption(self):
  c=certificate("none",1); self.assertTrue(verify(json.loads(json.dumps(c))))
  x=json.loads(json.dumps(c)); x["nodes"].append({"weights":[1]*8,"horizon":0,"value":"0","terminal":True,"actions":[]}); self.assertFalse(verify(x))
  x=json.loads(json.dumps(c)); n=next(n for n in x["nodes"] if n["actions"]); n["actions"]=n["actions"][:3]; self.assertFalse(verify(x))
  x=json.loads(json.dumps(c)); n=next(n for n in x["nodes"] if n["actions"]); n["actions"][0]["branches"][0]["probability"]="1/3"; self.assertFalse(verify(x))
  x=json.loads(json.dumps(c)); x["nodes"].append(json.loads(json.dumps(x["nodes"][0]))); self.assertFalse(verify(x))
  x=json.loads(json.dumps(c)); x["horizon"]="1"; self.assertFalse(verify(x))
