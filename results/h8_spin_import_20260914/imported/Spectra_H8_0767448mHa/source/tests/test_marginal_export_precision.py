import json,tempfile,unittest
from fractions import Fraction as F
from pathlib import Path
from experiments.marginal_hubbard_polynomial import build
from experiments.marginal_joint_coefficient_constructor import export,input_digest
from experiments.marginal_polynomial_metric import replay


class ExportPrecisionTests(unittest.TestCase):
    def test_exact_replay_at_bounded_higher_precision(self):
        raw=build(2,4,F(1));data={k:raw[k] for k in ('modes','particles','hamiltonian')}
        data['hamiltonian']=[t for t in data['hamiltonian'] if len(t['word'])==4]
        proposal={'source_sha256':input_digest(data),'gamma':'-14','orbits':[[[0,0]]],
            'metric_coefficients':[1.0000000001],'metric_bound':1.,'numerator_bound':18.,
            'weight_labels':[['charge',0,0]],'weight_values':[1.],
            'numerator_labels':[['charge',0,0]],'numerator_values':[18.]}
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder)
            for name,kwargs in [('default',{}),('fine',{'metric_precision':10**12,'proof_denominator':10**14})]:
                receipt=export(data,proposal,root/name,**kwargs)
                certificate=json.loads((root/name/'certificate.json').read_text())
                self.assertEqual(replay(certificate),receipt)
                self.assertGreater(F(receipt['weight_positivity']['lower']),0)
                self.assertGreater(F(receipt['numerator_positivity']['lower']),0)
            a=json.loads((root/'default/construction.json').read_text())
            b=json.loads((root/'fine/construction.json').read_text())
            self.assertEqual(F(a['exact_mean_Q_metric']),1)
            self.assertEqual(F(b['exact_mean_Q_metric']),F('1.0000000001'))
            for kwargs in ({'metric_precision':0},{'metric_precision':5*10**13+1},{'metric_precision':True},{'proof_denominator':10**16+1}):
                with self.assertRaises(ValueError):export(data,proposal,root/'invalid',**kwargs)


if __name__=='__main__':unittest.main()
