"""Bounded SU(2)-invariant Gram discovery with original-H exact replay."""
from fractions import Fraction as F
from pathlib import Path
from math import lcm,isfinite
import argparse,json,time
import numpy as np
import cvxpy as cp
from experiments.marginal_coefficient import dictionaries
from experiments.marginal_symbolic import canonical,hermitian,mono,add,scale,product,number_shift,multiplier_basis,encode,decode,verify
from experiments.marginal_hunt_car import adj
from research.certificate_scaling.spin_basis import decompose_words,spin_weight2
from research.certificate_scaling.spin_parity import spin_parity
from research.certificate_scaling.direct_sparse_discovery import sparse_columns
from research.certificate_scaling.spin_gram_columns import columns


def export(groups,grams,rounding=10**9):
    allrows=[];den=1;negative=0.
    if len(groups)!=len(grams):raise ValueError('Gram count mismatch')
    for group,Q in zip(groups,grams):
        k=len(group['copies'])
        if np.shape(Q)!=(k,k) or not np.all(np.isfinite(Q)):raise ValueError('Invalid Gram')
        vals,vecs=np.linalg.eigh((Q+Q.T)/2);negative+=float(-vals[vals<0].sum())
        raw=np.sqrt(np.maximum(vals,0))[:,None]*vecs.T*rounding
        if np.max(np.abs(raw),initial=0)>=2**62:raise ValueError('Excessive factor coefficient')
        Z=np.rint(raw).astype(np.int64)
        for r,weight in enumerate(group['weights']):
            rows=[]
            for row in Z:
                p=add(*(scale(copy[r],F(int(c),rounding)) for c,copy in zip(row,group['copies']) if c))
                if p:
                    rows.extend([p]*weight)
                    for c in p.values():den=lcm(den,c.denominator)
            if rows:allrows.append(rows)
    blocks=[]
    for rows in allrows:
        words=sorted({w for p in rows for w in p},key=lambda w:(len(w),w))
        blocks.append({'words':words,'factor':[[int(p.get(w,F(0))*den) for w in words] for p in rows]})
    return blocks,den,negative


def run(h,m,n,out,solver='CLARABEL',seconds=60):
    start=time.monotonic();out=Path(out)
    if type(m) is not int or m<=0 or m%2 or type(n) is not int or not 0<=n<=m:raise ValueError('Invalid paired-mode particle sector')
    if solver not in ('CLARABEL','SCS') or not isfinite(seconds) or seconds<=0:raise ValueError('Invalid solver budget')
    if not hermitian(h) or any(sum(2*c-1 for c,_ in w) for w in h):raise ValueError('Hermitian number-conserving H required')
    out.mkdir(parents=True,exist_ok=False)
    masks=spin_parity(h,m);groups=[];decompositions=[]
    for family in dictionaries(m,'mixed'):
        words=[w for w in family['words'] if len(w)!=1]
        if not words:continue
        gs,record=decompose_words(words,m,masks)
        for g in gs:g['family']=family['name']
        groups.extend(gs);decompositions.append({'family':family['name'],**record})
    built_basis=time.monotonic()
    def allowed(w):return spin_weight2(w)==0 and all(sum((mask>>i)&1 for c,i in w)%2==0 for mask in masks)
    basis=[p for p in multiplier_basis(m,2) if all(allowed(w) for w in p)]
    freepolys=[mono(())]+[product(number_shift(m,n),p) for p in basis]
    allwords=set(h)|{w for p in freepolys for w in p};polys=[];indices=[];work=0
    for group in groups:
        cs,ix,cost=columns(group);polys.append(cs);indices.append(ix);work+=cost
        allwords.update(w for p in cs for w in p)
    if any(len(w)>6 for w in allwords):raise AssertionError('Unexpected degree')
    rows=sorted(allwords,key=lambda w:(len(w),w));lookup={w:i for i,w in enumerate(rows)}
    free=sparse_columns(freepolys,lookup);maps=[sparse_columns(cs,lookup) for cs in polys]
    rhs=np.array([float(h.get(w,0)) for w in rows])
    rownorm=np.asarray(free.power(2).sum(axis=1)).ravel()
    for mat in maps:rownorm+=np.asarray(mat.power(2).sum(axis=1)).ravel()
    row_scale=1/np.maximum(1.,np.sqrt(rownorm));free_scale=1/np.maximum(1.,np.sqrt(np.asarray(free.power(2).sum(axis=0)).ravel()))
    x=cp.Variable(len(freepolys));rem=cp.Variable(len(rows));expr=free@cp.multiply(free_scale,x);qs=[]
    for g,mat,ix in zip(groups,maps,indices):
        k=len(g['copies']);q=cp.Variable((k,k),PSD=True);qs.append(q)
        expr+=mat@cp.reshape(q,(k*k,),order='C')[ix]
    problem=cp.Problem(cp.Minimize(cp.norm1(rem)-free_scale[0]*x[0]),[cp.multiply(row_scale,expr+rem-rhs)==0,rem[lookup[()]]==0])
    built=time.monotonic()
    metadata={'modes':m,'particles':n,'spin_parity_masks':masks,'decompositions':decompositions,'Gram_dimensions':[len(g['copies']) for g in groups],
              'Gram_entries':sum(len(g['copies'])**2 for g in groups),'coefficient_rows':len(rows),'gram_map_nonzeros':sum(mat.nnz for mat in maps),
              'polynomial_word_pair_products':work,'ideal_variables':len(basis),'basis_seconds':built_basis-start,'construction_seconds':built-start,
              'solver':solver,'solver_seconds_budget':seconds,'source_factors_used':False,'source_upper_used':False,
              'scope':'Invariant SOS in full chosen mixed-cubic word span. Original H retained; full exact residual is charged. No fixed-accuracy scaling theorem.'}
    (out/'pre_solve.json').write_text(json.dumps(metadata,indent=2)+'\n')
    options=({'tol_gap_abs':1e-8,'tol_feas':1e-8,'tol_gap_rel':1e-8,'max_iter':1000,'time_limit':float(seconds)} if solver=='CLARABEL' else {'eps':1e-8,'max_iters':150000,'time_limit_secs':float(seconds)})
    problem.solve(solver=solver,**options);solved=time.monotonic()
    if x.value is None or any(q.value is None for q in qs):raise RuntimeError('No numerical proposal')
    blocks,den,negative=export(groups,[q.value for q in qs]);xfree=x.value*free_scale
    X=add(*(scale(p,F(round(float(v)*10**12),10**12)) for p,v in zip(basis,xfree[1:])))
    cert={'modes':m,'particles':n,'hamiltonian':encode(h),'b':str(F(round(float(xfree[0])*10**12),10**12)),
          'number_multiplier':encode(X),'denominator':den,'blocks':blocks}
    exact=verify(cert);raw=json.dumps(cert,separators=(',',':'))+'\n';(out/'certificate.json').write_text(raw)
    result={**metadata,'status':problem.status,'solver_options':options,'solver_iterations':problem.solver_stats.num_iters,
            'solver_reported_seconds':problem.solver_stats.solve_time,'numeric_lower':-float(problem.value),'raw_b':float(xfree[0]),'negative_eigenvalue_mass':negative,
            'solve_seconds':solved-built,'export_seconds':time.monotonic()-solved,'wall_seconds':time.monotonic()-start,'exact':exact,'certificate_bytes':len(raw.encode())}
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ('decompositions','Gram_dimensions')}),flush=True)
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--fixture',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--solver',choices=['CLARABEL','SCS'],default='CLARABEL');p.add_argument('--seconds',type=float,default=60);a=p.parse_args();f=json.loads(a.fixture.read_text())
    run(decode(f['hamiltonian'],f['modes'],4),f['modes'],f['particles'],a.out,a.solver,a.seconds)
