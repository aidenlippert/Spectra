import json
from pathlib import Path
import tempfile
import unittest
import numpy as np
from scipy import sparse
from experiments.marginal_coefficient import dagger
from research.collective_completion_20260914.atoms import paired,initial,add,from_guide


class AtomTests(unittest.TestCase):
    def groups(self):
        words=[((1,0),(0,2),(0,1)),((1,2),(0,3),(0,0)),((1,4),(0,5),(0,0))]
        return [{'name':'cubic','words':words},{'name':'adjoint','words':[dagger(w) for w in words]}]
    def test_dual_violation_adds_a_positive_square(self):
        groups=self.groups();V=np.eye(3)[:,:1];members=paired(groups,[(0,V),(1,V)],V)
        maps=[(members,'initial')];meta={'pairs':[{'members':[0,1]}]};C=np.diag([1.,-2.,3.])
        with tempfile.TemporaryDirectory() as td:
            sparse.save_npz(Path(td)/'paired_map_0.npz',sparse.csc_matrix(C.reshape(1,-1)))
            result,ids,receipt=add(meta,groups,maps,[0],td,np.ones(1),1)
        self.assertEqual(ids,[0,0]);v=result[1][0][0][1]
        self.assertLess(float((v.T@C@v)[0,0]),0)
        np.testing.assert_array_equal(result[1][0][1][1],v)
    def test_seed_embedding_and_absent_zero_block(self):
        groups=self.groups();V=np.eye(3);maps=[(paired(groups,[(0,V),(1,V)],V),'cubic')]
        meta={'kind':'paired_quartic_v1','fixture_sha256':'same','rows':[],'pairs':[{'members':[0,1]}]}
        old={**meta,'groups':[{**g,'words':g['words'][:2]} for g in groups]}
        with tempfile.TemporaryDirectory() as td:
            p=Path(td);(p/'frame.json').write_text(json.dumps(old));(p/'construction.json').write_text(json.dumps({'prepared_dependency':td}))
            block={'name':'cubic','words':groups[0]['words'][:2],'factor':[[2,3]]}
            (p/'certificate.json').write_text(json.dumps({'core':{'denominator':1,'blocks':[block]}}))
            result,ids,receipt=initial(meta,groups,maps,p/'raw.npz')
            v=result[0][0][0][1].ravel();weight=receipt[0]['old_point_weights'][0]
            np.testing.assert_allclose(weight*np.outer(v,v),np.outer([2,3,0],[2,3,0]),atol=0)
            (p/'certificate.json').write_text(json.dumps({'core':{'denominator':1,'blocks':[]}}))
            result,ids,receipt=initial(meta,groups,maps,p/'raw.npz')
            self.assertEqual(receipt[0]['old_point_weights'],[0.]);self.assertEqual(len(result),1)
    def test_cold_guide_needs_no_certificate(self):
        groups=self.groups();V=np.eye(3)[:,:2];maps=[(paired(groups,[(0,V),(1,V)],V),'guide')]
        result,ids,receipt=from_guide(groups,maps)
        self.assertEqual(ids,[0,0]);self.assertEqual(receipt[0]['guide_atoms'],2)
        np.testing.assert_array_equal(np.hstack([m[0][0][1] for m in result]),V)


if __name__=='__main__':unittest.main()
