import contextlib, io, json, tempfile, unittest
from pathlib import Path
from fractions import Fraction as F
from experiments.marginal_valence_reference import build
from experiments.marginal_spin_temple import replay
from experiments.marginal_symbolic import mono,add,adj,product,scale,encode


def fixture():
    number=[mono(((1,i),(0,i))) for i in range(4)]
    h=add(scale(product(number[0],number[1]),4),scale(product(number[2],number[3]),4))
    for spin in (0,1):
        hop=mono(((1,spin),(0,2+spin)),F(-1,3));h=add(h,hop,adj(hop))
    return {'modes':4,'particles':2,'hamiltonian':encode(h)}


class ValenceReferenceTests(unittest.TestCase):
    def test_atom_gap_transfers_to_new_hopping_without_q_factor(self):
        from experiments.marginal_spin_transfer import transfer
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);source=root/'source.json';source.write_text(json.dumps(fixture()))
            with contextlib.redirect_stdout(io.StringIO()):
                build(source,root/'reference',upper_steps=2)
                transfer(root/'reference/certificate.json',root/'transfer',F(1,50),2,F(1,10**8),norm_squares=True,temple_offset=F(1,100))
            r=json.loads((root/'transfer/receipt.json').read_text());c=json.loads((root/'transfer/certificate.json').read_text())
            actual=replay(c)
            self.assertEqual(actual['width'],r['width']);self.assertLess(F(r['width']),F(1,10**8))
            self.assertEqual(r['largest_factor_dimension'],0)
            self.assertEqual(actual['complement_coverage']['reference_perturbation_norm_bound'],'1/25')
            self.assertEqual(actual['complement_coverage']['reference_coverage']['proof_family'],'spin_rational_atom_complement_v1')

    def test_full_small_interval_uses_empty_dd_gap_and_own_singlet_upper(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);source=root/'source.json';source.write_text(json.dumps(fixture()))
            with contextlib.redirect_stdout(io.StringIO()): result=build(source,root/'out',upper_steps=2)
            cert=json.loads((root/'out/certificate.json').read_text());r=replay(cert)
            self.assertLess(F(r['width']),F(1,10**8));self.assertEqual(result['retained_dimension'],2)
            self.assertEqual(result['gap_explicit_atom_count'],0);self.assertFalse(result['gap_has_dense_factor'])
            self.assertGreater(F(r['spin_symmetry_error_bound']),0)
            self.assertEqual(r['lower'],result['lower']);self.assertEqual(r['upper'],result['upper'])
            with self.assertRaises(ValueError): build(source,root/'out')
    def test_precision_and_budget_gates_precede_writing(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);source=root/'source.json';source.write_text(json.dumps(fixture()))
            for kwargs in ({'denominator':10},{'denominator':True},{'upper_steps':33},{'offset':F(0)}):
                with self.assertRaises(ValueError): build(source,root/'out',**kwargs)
                self.assertFalse((root/'out').exists())

if __name__=='__main__':unittest.main()
