"""Real inaccurate-iterate resume and canonical-identity refusal controls."""
import json,tempfile,unittest
from unittest.mock import patch
from pathlib import Path
import cvxpy as cp
from research.certificate_scaling.scs_checkpoint import solve,state_paths,identity


def problem(diagonal=1.):
    x=cp.Variable((2,2),PSD=True)
    return cp.Problem(cp.Minimize(-x[0,1]),[cp.diag(x)==diagonal]),x

class SCSCheckpointControls(unittest.TestCase):
    def test_backend_override_content_changes_identity(self):
        p,_=problem();data,_,_=p.get_problem_data(cp.SCS)
        with tempfile.TemporaryDirectory() as root:
            library=Path(root)/'experimental.so';library.write_bytes(b'first')
            with patch.dict('os.environ',{'SPECTRA_EVD_LIBRARY':str(library)}):
                first=identity(data,{},None);library.write_bytes(b'second')
                self.assertNotEqual(first,identity(data,{},None))
            with patch.dict('os.environ',{'SPECTRA_EVD_LIBRARY':'unresolved-library.so'}):
                with self.assertRaisesRegex(ValueError,'explicit library'):identity(data,{},None)

    def test_resume_inaccurate_iterate_and_change_budget(self):
        with tempfile.TemporaryDirectory() as root:
            root=Path(root);p,x=problem()
            first=solve(p,{'eps':1e-8,'max_iters':2},root/'first',contract={'case':'psd'})
            self.assertEqual(first['status_val'],2)
            # A fresh canonical problem, not a surviving CVXPY object/cache.
            q,y=problem()
            second=solve(q,{'eps':1e-8,'max_iters':2000},root/'second',root/'first',{'case':'psd'})
            self.assertEqual(second['resumed_from'],str(root/'first'))
            self.assertEqual(second['cumulative_iterations'],first['iterations_this_stage']+second['iterations_this_stage'])
            self.assertAlmostEqual(q.value,-1.,places=6)
            cold,z=problem();cold.solve(solver='SCS',eps=1e-8,max_iters=2000)
            self.assertAlmostEqual(q.value,cold.value,places=6)

    def test_changed_data_and_corruption_refused(self):
        with tempfile.TemporaryDirectory() as root:
            root=Path(root);p,x=problem();opts={'eps':1e-8,'max_iters':2}
            solve(p,opts,root/'first',contract={'case':'psd'})
            changed,_=problem(1.01)
            with self.assertRaisesRegex(ValueError,'identity'):solve(changed,opts,root/'bad',root/'first',{'case':'psd'})
            same,_=problem()
            with self.assertRaisesRegex(ValueError,'identity'):solve(same,opts,root/'bad',root/'first',{'case':'different'})
            meta,arrays=state_paths(root/'first');arrays.write_bytes(arrays.read_bytes()+b'corrupt')
            with self.assertRaisesRegex(ValueError,'hash'):solve(same,opts,root/'bad',root/'first',{'case':'psd'})

if __name__=='__main__':unittest.main()
