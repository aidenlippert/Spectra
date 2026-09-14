import copy
import json
from pathlib import Path
import tempfile
import unittest
from fractions import Fraction as F
from experiments.marginal_charge_metric_dual import propose,replay
from tests.test_marginal_charge_spin import fixed_sign_fixture


class ChargeMetricDualTests(unittest.TestCase):
    def test_physical_family_obstructions_and_false_target_refusal(self):
        with tempfile.TemporaryDirectory() as directory:
            p=Path(directory);source=p/'source.json';source.write_text(json.dumps(fixed_sign_fixture()))
            for family,classes in [('doublon_only',2),('one_pair_constant_tail',13)]:
                result=propose(source,p/family,family,'10');c=json.loads((p/family/'certificate.json').read_text())
                self.assertEqual(result,replay(c));self.assertEqual(result['metric_classes'],classes)
                self.assertTrue(all(F(x)<0 for x in result['class_coefficients']))
                bad=copy.deepcopy(c);bad['target_lower']='-100'
                with self.assertRaises(ValueError):replay(bad)
                for mutate in [lambda c:c['row_weights'][0].update(weight='-1'),
                               lambda c:c['row_weights'][0].update(state=0),
                               lambda c:c['row_weights'].append(c['row_weights'][0]),
                               lambda c:c.update(metric_family='unknown')]:
                    bad=copy.deepcopy(c);mutate(bad)
                    with self.assertRaises(ValueError):replay(bad)

    def test_no_obstruction_when_simple_positive_metric_is_feasible(self):
        with tempfile.TemporaryDirectory() as directory:
            p=Path(directory);source=p/'source.json';source.write_text(json.dumps(fixed_sign_fixture()))
            with self.assertRaises(ValueError):propose(source,p/'false','doublon_only','-100')


if __name__=='__main__':unittest.main()
