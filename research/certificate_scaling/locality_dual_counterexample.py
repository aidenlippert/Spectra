"""Exact CAR restricted-dual counterexample; full Fock sector, no number ideal.

A gapped vacuum does not constrain nonphysical duals of a monomial-square LP.
This refutes a gap-only pricing implication, not a theorem for stronger cones.
"""
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import argparse,json,sys,time
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from experiments.marginal_symbolic import word_product,add,mono,scale,canonical

def dagger(w):return tuple((1-c,i) for c,i in reversed(w))

def replay(m,t=F(1,4)):
    start=time.monotonic()
    if m<3 or not 0<t<F(1,2):raise ValueError('Gapped path requires M>=3, 0<t<1/2')
    chosen={(i,i+1) for i in range(m-1)}|{(i+1,i) for i in range(m-1)}|{(0,m-1),(m-1,0)}
    def y(poly):
        return sum(c for w,c in poly.items() if w==() or (len(w)==2 and w[0][0]==1 and w[1][0]==0 and (w[0][1],w[1][1]) in chosen))
    words=[tuple((1,i) for i in a)+tuple((0,i) for i in b) for d in range(3) for nc in range(d+1)
           for a in combinations(range(m),nc) for b in combinations(range(m),d-nc)]
    values=[y(dict(word_product(dagger(w),w))) for w in words]
    assert all(v>=0 for v in values)
    h=add(*(mono(((1,i),(0,i))) for i in range(m)),
          *(scale(add(mono(((1,i),(0,i+1))),mono(((1,i+1),(0,i)))),-t) for i in range(m-1)))
    hv=y(h);assert hv==-2*t*(m-1)
    endpoints=[((0,0),),((0,m-1),)]
    q=[[y(dict(word_product(dagger(u),v))) for v in endpoints] for u in endpoints]
    det=q[0][0]*q[1][1]-q[0][1]*q[1][0]
    ray=(q[0][0]-q[0][1]-q[1][0]+q[1][1])/F(2)
    assert q==[[0,1],[1,0]] and det==-1 and ray==-1
    # A primal certificate attaining the restricted optimum: sum_i n_i
    # is retained SOS; nearest-hop residual has l1 norm2t(M-1), b=0.
    residual=add(h,*(scale(mono(((1,i),(0,i))),-1) for i in range(m)))
    assert -sum(map(abs,residual.values()))==hv
    return {'modes':m,'graph_distance':m-1,'t':str(t),'retained_monomial_squares_checked':len(words),
            'retained_minimum':str(min(values)),'normalization':str(y({():F(1)})),
            'coefficient_dual_linf_bound':'1','restricted_optimum':str(hv),'primal_matches_dual':True,
            'endpoint_moment':[[str(x) for x in r] for r in q],'endpoint_determinant':str(det),
            'normalized_endpoint_rayleigh':str(ray),'analytic_full_Fock_gap_lower':str(1-2*t),
            'number_ideal_imposed':False,'far_atom_objective_gain_claimed':False,
            'wall_seconds':time.monotonic()-start}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--sizes',nargs='+',type=int,default=[4,8,12,16]);p.add_argument('--output',type=Path,default=Path('results/certificate_scaling/locality_dual_exact.json'));a=p.parse_args()
    rows=[replay(m) for m in a.sizes];a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(rows,indent=2)+'\n');print(json.dumps(rows,indent=2))
