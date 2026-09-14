import unittest
from fractions import Fraction as F
from copy import deepcopy
import json
from pathlib import Path
from experiments.v7_certificate import Generator,check_certificate,norm_witness,Piece
from experiments.v7_headroom import model,TOL
from experiments.v8_fractional_cover import check_fractional_cover,check_evolution_cover
from experiments.v11_guarded_taylor import guarded_taylor
from experiments.v12_small_cover import small_cover
from experiments.v12_small_cover_taylor import small_cover_taylor


class SmallCoverTests(unittest.TestCase):
    def test_fixed_residual_bound_and_exact_reconstruction(self):
        raw=json.loads((Path(__file__).resolve().parents[1]/'results/v8/fractional_probe.json').read_text())['rows'][10]
        op={p:F(c) for p,c in raw['residual'].items()}
        parts={mode:norm_witness(op,mode)[0] for mode in ('weighted','firstfit')}
        w,c=small_cover(op,parts)
        receipt=check_fractional_cover(op,w)
        self.assertEqual(receipt['status'],'certified')
        self.assertLessEqual(F(receipt['bound'])*F(raw['weight']),TOL)
        self.assertLessEqual(len(w['groups']),14)
        reconstructed={p:F(0) for p in op}
        for group in w['groups']:
            for p,value in group['coefficients'].items():
                weight=F(value)/op[p]
                self.assertGreaterEqual(weight,0);self.assertLessEqual(weight.denominator,1<<20)
                self.assertEqual((1<<20)%weight.denominator,0)
                reconstructed[p]+=F(value)
        self.assertEqual(reconstructed,op)
        bad=deepcopy(w);bad['groups'][0]['upper']='0'
        self.assertEqual(check_fractional_cover(op,bad)['status'],'rejected')
        bad=deepcopy(w);key=next(iter(bad['groups'][0]['coefficients']));bad['groups'][0]['coefficients'][key]=str(F(bad['groups'][0]['coefficients'][key])+1)
        self.assertEqual(check_fractional_cover(op,bad)['status'],'rejected')

    def test_complete_lower_degree_and_residual_sign(self):
        h,o=model(3,'xxz');T=F(1,2)
        p,w,c=small_cover_taylor(Generator(h,2,3,512),o,T,TOL)
        p0,w0,_=guarded_taylor(Generator(h,2,3,512),o,T,TOL)
        self.assertEqual(len(p.coefficients)-1,11);self.assertEqual(len(p0.coefficients)-1,12)
        self.assertEqual(p.coefficients,p0.coefficients[:12])
        ans=check_evolution_cover(Generator(h,2,3,512),o,[p],w,TOL,expected_time=T)
        self.assertEqual(ans['status'],'certified')
        bad=deepcopy(w)
        for group in bad['witnesses']['residual:0:11']['groups']:
            group['coefficients']={p:str(-F(v)) for p,v in group['coefficients'].items()}
        self.assertEqual(check_evolution_cover(Generator(h,2,3,512),o,[p],bad,TOL,expected_time=T)['status'],'rejected')

    def test_disabled_and_ordinary_success_match_baseline(self):
        h,o=model(3,'xxz');T=F(1,2)
        p,w,c=small_cover_taylor(Generator(h,2,3,512),o,T,TOL,enabled=False)
        p0,w0,_=guarded_taylor(Generator(h,2,3,512),o,T,TOL)
        self.assertEqual(p,p0);self.assertEqual(w['witnesses'],w0['witnesses']);self.assertEqual(c['cover_calls'],[])
        p,w,c=small_cover_taylor(Generator({},0,1),{'Z':F(1)},F(1),F(0))
        self.assertEqual(check_certificate(Generator({},0,1),{'Z':F(1)},[p],w,0,expected_time=1)['status'],'certified')

    def test_budget_and_partition_refusals(self):
        op={'X':F(1),'Z':F(1)}
        parts={m:norm_witness(op,m)[0] for m in ('weighted','firstfit')}
        with self.assertRaises(ValueError):small_cover(op,parts,grid_bits=33)
        bad=deepcopy(parts);bad['weighted'][0]['labels']=['X','X']
        with self.assertRaises(ValueError):small_cover(op,bad)
        w,_=small_cover({'X':F(0)},{'weighted':[],'firstfit':[]})
        self.assertEqual(check_fractional_cover({},w)['status'],'certified')


if __name__=='__main__':unittest.main()
