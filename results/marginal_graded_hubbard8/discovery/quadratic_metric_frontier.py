"""Reproducible finite 1106-pattern quadratic charge-metric frontier scan."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
import json
from itertools import product
from math import comb, prod
import numpy as np
from scipy.optimize import linprog
from experiments.marginal_joint_coefficient_constructor import feature_orbits

sites=8; orbits=feature_orbits(sites,2); cache={}
def f(q):
    if q not in cache:
        d=sum(x*x for x in q)//2
        cache[q]=np.array([d*sum(prod(x**p for x,p in zip(q,z)) for z in o) for o in orbits])
    return cache[q]
patterns=[q for q in product((-1,0,1),repeat=sites) if sum(q)==0 and any(q)]
W=[];K=[];mult=[]
for q in patterns:
    v=f(q);d=sum(x*x for x in q)//2;pen=np.zeros(len(orbits))
    for i in range(sites-1):
        a,b=q[i]+1,q[i+1]+1
        if abs(a-b)==1:
            t=list(q);t[i],t[i+1]=t[i+1],t[i];pen+=f(tuple(t))
        elif {a,b}=={0,2}:
            t=list(q);t[i]=t[i+1]=0;pen+=2*f(tuple(t))
        elif a==b==1:
            for s in (-1,1):
                t=list(q);t[i],t[i+1]=s,-s;pen+=f(tuple(t))
    W.append(v);K.append(4*d*v-pen);mult.append(comb(sites-2*d,(sites-2*d)//2))
W,K=np.array(W),np.array(K); mean=np.array(mult)@W/sum(mult)
rows=[]
for gamma in (-12,-10,-8,-6,-4):
    # maximize t subject to Kc-gamma Wc >= t and Wc >= 0, mean*c=1
    A=np.vstack((np.c_[-(K-gamma*W),np.ones(len(patterns))],np.c_[-W,np.zeros(len(patterns))]))
    r=linprog(np.r_[np.zeros(len(orbits)),-1],A_ub=A,b_ub=np.zeros(2*len(patterns)),A_eq=[np.r_[mean,0]],b_eq=[1],bounds=[(None,None)]*(len(orbits)+1),method='highs-ds')
    row={'gamma':gamma,'success':bool(r.success),'maximum_numerator_floor':float(r.x[-1]) if r.success else None,'minimum_metric':float(min(W@r.x[:-1])) if r.success else None,'status':r.message}
    if gamma == -10 and r.success:
        # Preserve the raw numerical primal witness for a separate exact-dual
        # reconstruction attempt; it is intentionally not called a certificate.
        row['metric_coefficients']=r.x[:-1].tolist()
    rows.append(row)
out={'sites':8,'charge_patterns':len(patterns),'spin_sector_multiplicity_sum':sum(mult),'results':rows,'scope':'Finite numerical diagnostic over all 1106 neutral charge patterns. Nonnegative normalized quadratic charge metric family; negative floors are not exact infeasibility certificates.'}
dest=Path('results/marginal_graded_hubbard8/quadratic_metric_frontier');dest.mkdir(parents=True,exist_ok=True)
(dest/'receipt.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
