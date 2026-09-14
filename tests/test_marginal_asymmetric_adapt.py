from fractions import Fraction as F
import unittest
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import numpy as np
from experiments.marginal_asymmetric_adapt import exchange_subspaces,exchange_projection,prepare
from experiments.marginal_dual_exact import full_word_blocks
from experiments.marginal_reynolds import invariant_rows
from experiments.marginal_orbit_quartic import word_action,symmetry_group
from experiments.marginal_orbit_certificate import average_canonical_polynomial,verify
from experiments.marginal_symbolic import adj,add,mono,canonical,encode,product
from experiments.marginal_collective import hopping_polynomial


class AsymmetricAdaptTests(unittest.TestCase):
    def test_orbit_columns_and_certified_atom_match_source_algebra(self):
        p={((0,0),):F(1),((0,2),):F(2)}
        base={'modes':4,'particles':2,'operator_degree':4,'hamiltonian':encode(hopping_polynomial(4,F(1,5))),
              'b':'0','number_multiplier':[],'permutations':[list(g) for g in symmetry_group(4)],
              'orbit_squares':[{'polynomial':encode(p),'weight':'1'}]}
        with TemporaryDirectory() as directory:
            path=Path(directory)/'source.json';path.write_text(json.dumps(base))
            data=prepare(path);seed=data['seeds'][0]
            square=average_canonical_polynomial(product(adj(seed),seed),base['permutations'])
            expected=data['projection']@np.array([float(square.get(w,0)) for w in data['rows']])
            self.assertTrue(np.allclose(data['columns'][:,0],expected,atol=1e-13,rtol=0))
            atom=prepare(path,aggregate=True);lower=F(verify(base)['lower']);h=hopping_polynomial(4,F(1,5))
            expected=atom['projection']@np.array([float(h.get(w,0)-(lower if w==() else 0)) for w in atom['rows']])
            self.assertTrue(np.allclose(atom['columns'][:,0],expected,atol=1e-13,rtol=0))

    def test_exchange_subspaces_are_complete_eigenspaces(self):
        modes=6;mapping=[(i+3)%6 for i in range(6)]
        for block in full_word_blocks(modes,3):
            words=block['words'];lookup={w:i for i,w in enumerate(words)};operator=np.zeros((len(words),len(words)))
            for j,w in enumerate(words):
                image,sign=word_action(w,mapping);operator[lookup[image],j]=float(sign)
            subspaces=exchange_subspaces(words,modes);u=np.column_stack(subspaces)
            self.assertTrue(np.allclose(u@u.T,np.eye(len(words)),atol=1e-13,rtol=0))
            for tr in subspaces:
                self.assertTrue(np.allclose(operator@tr,tr,atol=1e-13,rtol=0) or np.allclose(operator@tr,-tr,atol=1e-13,rtol=0))

    def test_exchange_projection_matches_exact_average_and_residual_norm(self):
        modes=6;rows=invariant_rows(modes,4);projection,weights=exchange_projection(rows,modes)
        p={}
        for j,w in enumerate(rows):
            term=mono(w,F(j%7-3,11));p=add(p,term,adj(term))
        p=canonical(p)
        average=average_canonical_polynomial(p,[list(range(modes)),[(i+3)%6 for i in range(6)]])
        actual=projection@np.array([float(p.get(w,0)) for w in rows])
        expected=projection@np.array([float(average.get(w,0)) for w in rows])
        self.assertTrue(np.allclose(actual,expected,atol=1e-13,rtol=0))
        self.assertAlmostEqual(float(weights@abs(actual)),float(sum(abs(c) for c in average.values())),places=10)


if __name__=='__main__':unittest.main()
