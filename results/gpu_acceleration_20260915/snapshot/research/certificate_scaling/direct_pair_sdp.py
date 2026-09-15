"""H-local pair-block SDSOS proposal with a full CAR residual objective."""
from __future__ import annotations
import argparse, json, sys, time
from fractions import Fraction as F
from pathlib import Path
import numpy as np
import cvxpy as cp
from scipy import sparse
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from experiments.marginal_symbolic import add, canonical, decode, encode, multiplier_basis, mono, number_shift, product, scale, verify, word_product
from research.certificate_scaling.direct_sparse_discovery import candidates
from experiments.marginal_coefficient import dagger

def run(h, modes, particles, budget, outdir, denominator=10**6):
    started=time.monotonic(); h=canonical(h); atoms, counts=candidates(h,modes,budget)
    singles=[]; pairs=[]
    for ws, signs, _ in atoms:
        if len(ws)==1: singles.append(ws[0])
        elif len(ws)==2:
            key=tuple(sorted(ws));
            if key not in [p[0] for p in pairs]: pairs.append((key,))
    pairs=[p[0] for p in pairs]
    basis=multiplier_basis(modes,max_body=1); shift=number_shift(modes,particles); ideal=[product(shift,p) for p in basis]
    polys=[mono(())]+ideal; labels=[('b',None)]+[('x',i) for i in range(len(ideal))]
    for w in singles: polys.append(dict(word_product(dagger(w),w))); labels.append(('s',w))
    for u,v in pairs:
        polys += [dict(word_product(dagger(u),u)), dict(word_product(dagger(u),v)), dict(word_product(dagger(v),u)), dict(word_product(dagger(v),v))]
    words=sorted(set(h)|{w for p in polys for w in p},key=lambda w:(len(w),w)); lookup={w:i for i,w in enumerate(words)}
    def col(p): return np.array([float(p.get(w,0)) for w in words])
    b=cp.Variable(); x=cp.Variable(len(ideal)); rem=cp.Variable(len(words)); qsingle=cp.Variable(len(singles),nonneg=True); qpair=[cp.Variable((2,2),PSD=True) for _ in pairs]
    expr=b*col(polys[0])+sum(x[i]*col(ideal[i]) for i in range(len(ideal)))
    for i in range(len(singles)): expr += qsingle[i]*col(polys[1+len(ideal)+i])
    off=1+len(ideal)+len(singles)
    for i in range(len(pairs)):
        base=off+4*i; q=qpair[i]
        expr += q[0,0]*col(polys[base])+q[0,1]*col(polys[base+1])+q[1,0]*col(polys[base+2])+q[1,1]*col(polys[base+3])
    rhs=np.array([float(h.get(w,0)) for w in words]); constraints=[expr+rem==rhs]
    problem=cp.Problem(cp.Minimize(cp.norm1(rem)-b),constraints)
    build=time.monotonic()-started; solve0=time.monotonic(); problem.solve(solver='SCS',eps=1e-7,max_iters=5000,time_limit_secs=60); solve=time.monotonic()-solve0
    if b.value is None: raise RuntimeError(f'SDP failed: {problem.status}')
    blocks=[]
    for i,w in enumerate(singles):
        z=int(round(np.sqrt(max(float(qsingle.value[i]),0))*denominator));
        if z: blocks.append({'name':f'single-{i}','words':[w],'factor':[[z]]})
    for i,(u,v) in enumerate(pairs):
        vals,vec=np.linalg.eigh((qpair[i].value+qpair[i].value.T)/2); L=(np.sqrt(np.maximum(vals,0))[:,None]*vec.T); ints=np.rint(L*denominator).astype(int).tolist(); ints=[r for r in ints if any(r)]
        if ints: blocks.append({'name':f'pair-{i}','words':[list(u),list(v)],'factor':ints})
    xv=add(*(scale(p,F(int(round(float(x.value[i])*10**9)),10**9)) for i,p in enumerate(basis) if x.value is not None))
    cert={'modes':modes,'particles':particles,'hamiltonian':encode(h),'b':str(F(int(round(float(b.value)*10**9)),10**9)),'number_multiplier':encode(xv),'denominator':denominator,'blocks':blocks}
    receipt=verify(cert); done=time.monotonic(); receipt.update(counts|{'method':'direct_H_pair_SDP_l1','budget':budget,'pair_blocks':len(pairs),'singleton_blocks':len(singles),'coefficient_rows':len(words),'max_degree':max(map(len,words)),'gram_scalar_variables':len(singles)+4*len(pairs),'build_seconds':build,'solve_seconds':solve,'wall_seconds':done-started,'certificate_bytes':len(json.dumps(cert,separators=(',',':')).encode()),'source_factors_used':False,'source_upper_used':False,'sdp_status':problem.status})
    outdir.mkdir(parents=True,exist_ok=True); (outdir/'certificate.json').write_text(json.dumps(cert,separators=(',',':'))+'\n'); (outdir/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n'); return receipt

def main():
    p=argparse.ArgumentParser(); p.add_argument('--fixture',type=Path,required=True); p.add_argument('--budget',type=int,default=64); p.add_argument('--outputdir',type=Path,required=True); p.add_argument('--denominator',type=int,default=10**6); a=p.parse_args(); d=json.loads(a.fixture.read_text()); h=decode(d['hamiltonian'],d['modes'],4); print(json.dumps(run(h,d['modes'],d['particles'],a.budget,a.outputdir,a.denominator),indent=2))
if __name__=='__main__': main()
