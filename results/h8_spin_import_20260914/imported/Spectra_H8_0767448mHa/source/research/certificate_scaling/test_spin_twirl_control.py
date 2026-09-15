import contextlib,io,json,tempfile,unittest
from fractions import Fraction as F
from pathlib import Path
from experiments.marginal_symbolic import decode,expand_squares,add,scale,number_shift,product,mono
from research.certificate_scaling.spin_twirl import twirl
from research.certificate_scaling.spin_twirl_control import run

class TwirlControlTests(unittest.TestCase):
    def test_original_h4_control_matches_exact_projected_source_identity(self):
        source=Path(__file__).resolve().parents[2]/'results/certificate_scaling/cubic_precision/local_h4_quotient/certificate.json'
        old=json.loads(source.read_text());m=old['modes']
        old_square,_=expand_squares(old['blocks'],old['denominator'],m)
        with tempfile.TemporaryDirectory() as d,contextlib.redirect_stdout(io.StringIO()):
            result=run(source,Path(d)/'out');new=json.loads((Path(d)/'out/certificate.json').read_text())
        new_square,_=expand_squares(new['blocks'],new['denominator'],m)
        # The expected projection uses polynomial Casimir actions, independently
        # of local coordinate inversion and numerical Gram contraction.
        diff=add(new_square,scale(twirl(old_square),-1))
        self.assertLess(sum(abs(c) for c in diff.values()),F(1,1000000))
        self.assertEqual(new['hamiltonian'],old['hamiltonian'])
        self.assertTrue(result['source_certificate_used'])
        self.assertEqual(result['Gram_entries'],5088)
        self.assertLess(result['negative_eigenvalue_mass'],1e-10)
        self.assertGreater(F(result['exact']['lower']),F(-3667001,1000000))

if __name__=='__main__':unittest.main()
