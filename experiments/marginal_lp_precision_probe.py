"""Audit HiGHS-to-exact export loss for the asymmetric polishing LP.

This is deliberately a standalone probe; it does not modify the production
polisher or any existing result directory.
"""
from fractions import Fraction as F
import json, time
from pathlib import Path
import numpy as np
import scipy.sparse as sp
from scipy.optimize import linprog
from experiments.marginal_asymmetric_adapt import prepare
from experiments.marginal_reynolds import reynolds_matrix
from experiments.marginal_symbolic import decode, encode, product, adj, add, scale
from experiments.marginal_orbit_certificate import average_canonical_polynomial, normalize_scalar_residual
from experiments.marginal_symmetry_transfer import replay

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/marginal_lp_precision_probe'

def run(source):
    t0=time.monotonic(); source=Path(source); original=json.loads(source.read_text())
    data=prepare(ROOT/'results/marginal_distill_direct/compact_certificate.json')
    m=data['modes']; rows=data['rows']; projection=data['projection']; full_average=reynolds_matrix(rows,m)[0]
    seeds=[]; cols=[]; dropped=[]
    for kind in ('orbit_squares','direct_squares'):
        for item in original.get(kind,[]):
            p0=decode(item['polynomial'],m,4); largest=max(abs(c) for c in p0.values()); p=scale(p0,1/largest)
            sq=product(adj(p),p); col=np.array([float(sq.get(w,0)) for w in rows])
            if kind=='orbit_squares': col=full_average@col
            seeds.append((kind,p)); cols.append(projection@col)
    n=projection.shape[0]; rank=len(data['basis']); nr=len(seeds)
    matrix=sp.hstack([sp.csr_matrix(data['unit'][:,None]),sp.csr_matrix(data['a']),sp.csr_matrix(np.column_stack(cols)),sp.eye(n),-sp.eye(n)],format='csr')
    h=decode(original['hamiltonian'],m,4); rhs=projection@np.array([float(h.get(w,0)) for w in rows])
    objective=np.r_[-1,np.zeros(rank+nr),data['weights'],data['weights']]
    rb=[(0,0)]+[(0,None)]*(n-1); bounds=[(None,None)]*(1+rank)+[(0,100)]*nr+rb+rb
    t1=time.monotonic(); result=linprog(objective,A_eq=matrix,b_eq=rhs,bounds=bounds,method='highs',options={'primal_feasibility_tolerance':1e-9,'dual_feasibility_tolerance':1e-9}); lp=time.monotonic()-t1
    tag=source.parent.name+'_'+source.stem
    out=OUT/tag
    status={'source':str(source),'success':bool(result.success),'status':result.message,'build_seconds':t1-t0,'lp_seconds':lp,'rows':n,'candidate_squares':nr,'multiplier_rank':rank,'options_used':{'method':'highs','primal_feasibility_tolerance':1e-9,'dual_feasibility_tolerance':1e-9}}
    if not result.success:
        OUT.mkdir(exist_ok=True); (OUT/'status.json').write_text(json.dumps(status,indent=2)+'\n'); return status
    x=result.x; raw_res=np.asarray(matrix@x-rhs); status.update(numerical_objective=float(-result.fun),max_equation_residual=float(np.max(np.abs(raw_res))),equation_l2=float(np.linalg.norm(raw_res)),weighted_equation_l1=float(np.dot(data['weights'],np.abs(raw_res))),weighted_equation_linf=float(np.max(data['weights']*np.abs(raw_res))))
    nz=np.abs(matrix.data); status.update(matrix_nonzero_entries=int(nz.size),matrix_min_abs_nonzero=float(nz.min()),matrix_quantiles_abs_nonzero=[float(v) for v in np.quantile(nz,[0,.01,.5,.99,1])])
    def rat(c,d=10**12): return F(round(float(c)*d),d)
    selected={'orbit_squares':[],'direct_squares':[]}; kept_l1=0.; dropped_l1=0.; kept=0; drop=0
    for (kind,p), weight in zip(seeds,x[1+rank:1+rank+nr]):
        w=rat(max(0,weight)); rounded={word:F(round(c*10**14),10**14) for word,c in p.items() if round(c*10**14)}
        dropped_l1 += sum(abs(float(c)) for word,c in p.items() if not round(c*10**14)); kept_l1 += sum(abs(float(c)) for word,c in p.items() if round(c*10**14)); drop += sum(1 for c in p.values() if not round(c*10**14)); kept += len(rounded)
        if w and rounded: selected[kind].append({'polynomial':encode(rounded),'weight':str(w)})
    multiplier=average_canonical_polynomial(add(*(scale(p,rat(c)) for p,c in zip(data['basis'],x[1:1+rank]))),[list(range(m)),[(i+5)%m for i in range(m)]])
    cert={'modes':m,'particles':original['particles'],'operator_degree':4,'hamiltonian':original['hamiltonian'],'b':str(rat(x[0])),'number_multiplier':encode(multiplier),'permutations':original['permutations'],'independent_upper':original['independent_upper'],**selected}
    cert,_=normalize_scalar_residual(cert); receipt=replay(cert)
    status.update(exported_b=cert['b'],lower_float=receipt.get('lower_float'),width_float=receipt.get('width_float'),export_gap=float(status['numerical_objective']-receipt.get('lower_float',0)),polynomial_terms_kept=kept,polynomial_terms_dropped=drop,coefficient_l1_kept=kept_l1,coefficient_l1_dropped=dropped_l1,residual_l1_float=receipt.get('residual_l1_float'),residual_nonconstant_l1=receipt.get('residual_nonconstant_l1'),orbit_squares=receipt.get('orbit_squares'),direct_squares=receipt.get('direct_squares'),total_seconds=time.monotonic()-t0)
    out.mkdir(parents=True,exist_ok=True); (out/'status.json').write_text(json.dumps(status,indent=2)+'\n'); (out/'certificate.json').write_text(json.dumps(cert)+'\n'); (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n'); return status

if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser(); ap.add_argument('source',type=Path); print(json.dumps(run(ap.parse_args().source),indent=2))
