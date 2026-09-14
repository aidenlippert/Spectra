#!/usr/bin/env python3
"""Fit finite-range positive product metrics to the exact charge PF vector.

The emitted tables are rational positive factors (common denominator 1e8).  All
row bounds are evaluated from the jointly assembled product metric, so shared
factors are retained across transitions in one physical row.
"""
from itertools import product
from pathlib import Path
import json, math, hashlib
import numpy as np
from scipy.optimize import least_squares, minimize

ROOT=Path(__file__).resolve().parents[3]
CERT=ROOT/'results/marginal_graded_hubbard8/charge_only_limit/certificate.json'
OUT=ROOT/'results/marginal_graded_hubbard8/product_metric_fit'
GEN=ROOT/'results/marginal_graded_hubbard8/discovery/charge_only_limit.py'

def build(sites=8):
    ps=[q for q in product((-1,0,1),repeat=sites) if sum(q)==0 and any(q)]; ix={q:i for i,q in enumerate(ps)}; A={}
    for r,q in enumerate(ps):
        A[(r,r)]=-4*(sum(x*x for x in q)//2)
        for i in range(sites-1):
            a,b=q[i]+1,q[i+1]+1; zs=[]
            if abs(a-b)==1:
                z=list(q); z[i],z[i+1]=z[i+1],z[i]; zs.append((tuple(z),1))
            elif {a,b}=={0,2}:
                z=list(q); z[i]=z[i+1]=0; zs.append((tuple(z),2))
            elif a==b==1:
                for d in (-1,1):
                    z=list(q); z[i],z[i+1]=d,-d; zs.append((tuple(z),1))
            for z,w in zs:
                if z in ix:A[(r,ix[z])]=A.get((r,ix[z]),0)+w
    return ps,A

def fit(R):
    cert=json.loads(CERT.read_text()); ps,A=build(); w=np.array(cert['right_weight'],float); sites=8
    # intercept omitted; gauge fixes every onsite(0) and pair(0,0) to one.
    keys=[('o',i,c) for i in range(sites) for c in (-1,1)]
    keys += [('p',i,j,a,b) for i in range(sites) for j in range(i+1,min(sites,i+R+1)) for a in (-1,0,1) for b in (-1,0,1) if (a,b)!=(0,0)]
    X=np.zeros((len(ps),len(keys)))
    for r,q in enumerate(ps):
        for k,key in enumerate(keys):
            if key[0]=='o': X[r,k]=1 if q[key[1]]==key[2] else 0
            else: X[r,k]=1 if (q[key[1]],q[key[2]])==(key[3],key[4]) else 0
    D=np.array([sum(x*x for x in q)//2 for q in ps],float)
    # Overall metric scale is immaterial; remove the global log offset to make
    # the gauge-fixed factorization numerically well conditioned.
    y=np.log(w/D); y=y-y.mean()
    ls=least_squares(lambda z:X@z-y,np.zeros(len(keys)),bounds=(-8,8),max_nfev=3000)
    z=ls.x
    # Direct smooth-minimax refinement of the physical row bounds.  For each
    # row g=4D-sum rate exp((X_c-X_r)z), with analytic gradient.
    trans=[[(c,rate) for (rr,c),rate in A.items() if rr==r and c!=r] for r in range(len(ps))]
    def gg(x):
        g=4*D.copy(); J=np.zeros((len(ps),len(keys)))
        for r,ts in enumerate(trans):
            for c,rate in ts:
                rat=rate*np.exp(np.clip((X[c]-X[r])@x,-60,60)); g[r]-=rat; J[r]-=rat*(X[c]-X[r])
            g[r]+=4*D[r]
        return g,J
    for tau in ():
        def fun(x):
            g,J=gg(x); a=-tau*g; m=a.max(); p=np.exp(a-m); p/=p.sum()
            return float(-(m+np.log(np.exp(a-m).sum()))/tau), -p@J
        opt=minimize(fun,z,jac=True,method='L-BFGS-B',bounds=[(-5,5)]*len(keys),options={'maxiter':1200,'ftol':1e-10})
        z=opt.x
    # A bounded least-squares fit is stable for this gauge-fixed, sparse basis.
    # (The exact row bound below is the acceptance objective.)
    """
    def obj(x):
        e=X@x-y; tau=40.; m=np.max(e); return m+np.log(np.exp(tau*(e-m)).sum())/tau
    opt=minimize(obj,z,method='L-BFGS-B',bounds=[(-8,8)]*len(keys),options={'maxiter':3000,'ftol':1e-12})
    z=opt.x
    """
    vals=np.exp(z); den=100_000_000; nums=np.maximum(1,np.rint(vals*den).astype(np.int64)); zq=np.log(nums/den)
    logv=np.log(D)+X@zq; v=np.exp(logv)
    rows=[]
    for r,q in enumerate(ps):
        out=sum(rate*v[c]/v[r] for (rr,c),rate in A.items() if rr==r and c!=r)
        rows.append(4*D[r]-out)
    payload={'range':R,'sites':sites,'patterns':len(ps),'source_generator':str(GEN.relative_to(ROOT)),'source_sha256':hashlib.sha256(GEN.read_bytes()).hexdigest(),'denominator':den,'factor_parameter_count':len(keys),'onsite':[[int(nums[keys.index(('o',i,c))]) if c else den for c in (-1,0,1)] for i in range(sites)],'pair':{},'fit_rmse':float(np.sqrt(np.mean((X@z-y)**2))),'fit_max_abs_log_error':float(np.max(np.abs(X@z-y))),'weighted_row_bound_min':float(min(rows)),'weighted_row_bound_max':float(max(rows)),'target_gamma':-4.0,'target_gamma_reached':bool(min(rows)>=-4.0),'certificate_pattern_digest':hashlib.sha256(CERT.read_bytes()).hexdigest()}
    for i in range(sites):
        for j in range(i+1,min(sites,i+R+1)):
            payload['pair'][f'{i},{j}']=[[int(nums[keys.index(('p',i,j,a,b))]) if (a,b)!=(0,0) else den for b in (-1,0,1)] for a in (-1,0,1)]
    payload['row_bound_units']='gamma=4D(q)-sum_transition_rate*v(target)/v(source)'; payload['best_integer_fraction_target']=float(math.floor(min(rows)))
    return payload

def main():
    OUT.mkdir(parents=True,exist_ok=True); allr={}
    for R in (1,2,3):
        p=fit(R); allr[str(R)]=p; (OUT/f'R{R}.json').write_text(json.dumps(p,indent=2)+'\n')
    (OUT/'summary.json').write_text(json.dumps({'status':'complete','ranges':allr,'scope':'exact charge-only matrix; no enumeration rediscovery claim'},indent=2)+'\n')
if __name__=='__main__': main()
