from fractions import Fraction as F
from itertools import combinations
import unittest
import numpy as np
from experiments.marginal_symbolic import mono,add,encode,canonical,adj,product,scale
from research.certificate_scaling.spin_twirl import twirl
from research.collective_completion_20260914.spin_rows import build
from research.collective_completion_20260914.spin_screen import check_sector

class SpinRowsTests(unittest.TestCase):
    def test_exact_projector_basis(self):
        rows=[tuple((1,i) for i in a)+tuple((0,i) for i in b) for k in range(3)
              for a in combinations(range(4),k) for b in combinations(range(4),k)
              if a<=b and sum(i%2==0 for i in a)==sum(i%2==0 for i in b)]
        T,selected,rec=build(rows)
        self.assertLess(len(selected),len(rows));self.assertTrue(rec['rank_equals_exact_trace'])
        coeff=np.arange(len(rows))-3;poly={}
        for w,c in zip(rows,coeff):
            a=canonical(adj(mono(w)));p=mono(w) if a==mono(w) else add(mono(w),a)
            poly=add(poly,scale(p,F(int(c))))
        exact=twirl(poly)
        np.testing.assert_allclose(T@coeff,[float(exact.get(w,0)) for w in rows],atol=1e-12)
    def test_averaged_positive_square(self):
        h=twirl(mono(((1,0),(0,0))));data={'modes':4,'particles':2,'hamiltonian':encode(h)}
        core={**data,'operator_degree':3,'number_multiplier':[],'b':'0','denominator':1,
              'blocks':[{'words':[[[0,0]]],'factor':[[1]]}]}
        cert={**data,'kind':'spin_sector_sos_v1','core':core,'alpha_multiplier':[],'casimir_multiplier':'0',
              'magnetization':0,'singlet':True,'spin_twirl':True}
        r=check_sector(data,cert);self.assertEqual(F(r['lower']),0);self.assertEqual(F(r['unaveraged_residual_l1']),1)
        bad=add(h,mono((),F(-1,10)));data['hamiltonian']=encode(bad);cert['hamiltonian']=data['hamiltonian'];core['hamiltonian']=data['hamiltonian']
        self.assertEqual(F(check_sector(data,cert)['lower']),F(-1,10))
        cert['singlet']=False
        with self.assertRaises(ValueError):check_sector(data,cert)
if __name__=='__main__':unittest.main()
