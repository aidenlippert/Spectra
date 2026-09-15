import json, tempfile, unittest
from pathlib import Path
from experiments.marginal_localized_search import construct, _dictionary
from experiments.marginal_signed_atoms import replay, atom_vector
from experiments.marginal_clique_gap import balanced_signs
from fractions import Fraction as F
from tests.test_marginal_localized_basis import fixture

class LocalizedSearchTests(unittest.TestCase):
    def test_empty_seed_replays_and_dictionary_is_bounded(self):
        data,_=fixture()
        with tempfile.TemporaryDirectory() as z:
            p=Path(z); src=p/'h.json'; out=p/'out'; src.write_text(json.dumps(data))
            # fixture has four spin-sector states; retain one physical determinant
            r=construct(src,out,[3],target='10')
            self.assertEqual(r['blocks'],1); self.assertFalse(r['requested_target_certified'])
            cert=json.loads((out/'certificate.json').read_text()); metric=json.loads((out/'metric_source.json').read_text())
            state=json.loads((out/'search_state.json').read_text())
            self.assertEqual(metric['target_lower'],'10'); self.assertEqual(state['vector_mode'],'rational')
            self.assertEqual(replay(cert)['target_lower'],r['target_lower'])
            self.assertLess(F(cert['target_lower']),F(metric['target_lower']))

    def test_dictionary_contains_only_balanced_local_supports(self):
        a=[[F(1) for j in range(9)] for i in range(9)]
        a[0][1]=a[1][0]=F(-1)
        entries=_dictionary(a)
        self.assertGreater(len(entries),0); self.assertLessEqual(len(entries),9*35)
        for item in entries:
            group,signs=atom_vector(item,9,True)
            self.assertEqual(list(signs),balanced_signs(a,group))
    def test_no_overwrite(self):
        data,_=fixture()
        with tempfile.TemporaryDirectory() as z:
            p=Path(z); src=p/'h.json'; out=p/'out'; src.write_text(json.dumps(data)); out.mkdir()
            with self.assertRaises(ValueError): construct(src,out,[3])

if __name__=='__main__': unittest.main()
