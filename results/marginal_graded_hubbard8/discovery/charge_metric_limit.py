"""Explicit small charge-pattern diagnostic, separate from implicit construction."""
from pathlib import Path
from itertools import product
from math import comb, prod
import json
import numpy as np
from scipy.optimize import linprog
from experiments.marginal_joint_coefficient_constructor import feature_orbits

sites = 8; orbits = feature_orbits(sites, 2); cache = {}
def features(q):
    if q not in cache:
        D = sum(x*x for x in q)//2
        cache[q] = np.array([D*sum(prod(x**p for x,p in zip(q,powers)) for powers in orbit) for orbit in orbits], dtype=float)
    return cache[q]
patterns = [q for q in product((-1,0,1), repeat=sites) if sum(q)==0 and any(q)]
weights=[]; numerator=[]; multiplicities=[]
for q in patterns:
    v=features(q); D=sum(x*x for x in q)//2; penalty=np.zeros(len(orbits))
    for i in range(sites-1):
        a,b=q[i]+1,q[i+1]+1
        if abs(a-b)==1:
            target=list(q);target[i],target[i+1]=target[i+1],target[i];penalty+=features(tuple(target))
        elif {a,b}=={0,2}:
            target=list(q);target[i]=target[i+1]=0;penalty+=2*features(tuple(target))
        elif a==b==1:
            for direction in (-1,1):
                target=list(q);target[i]=direction;target[i+1]=-direction;penalty+=features(tuple(target))
    weights.append(v);numerator.append(4*D*v-penalty);multiplicities.append(comb(sites-2*D,(sites-2*D)//2))
W,K=np.array(weights),np.array(numerator);mean=np.array(multiplicities)@W/sum(multiplicities);n=len(orbits)
receipts=[]
for gamma in (-18,-16,-14):
    A=np.vstack((np.c_[-(K-gamma*W),np.ones(len(patterns))],np.c_[-W,np.zeros(len(patterns))]))
    result=linprog(np.r_[np.zeros(n),-1.],A_ub=A,b_ub=np.r_[np.zeros(len(patterns)),np.zeros(len(patterns))],A_eq=np.array([np.r_[mean,0.]]),b_eq=[1.],bounds=[(None,None)]*(n+1),method='highs-ds')
    receipt={'gamma':gamma,'status':result.message,'success':bool(result.success)}
    if result.success:
        receipt.update(maximum_numerator_floor=float(result.x[-1]),minimum_metric=float(min(W@result.x[:n])),metric_coefficients=result.x[:n].tolist(),max_inequality_violation=float(max(A@result.x)))
    receipts.append(receipt)
out={'sites':sites,'charge_patterns':len(patterns),'spin_sector_multiplicity_sum':sum(multiplicities),'feature_orbits':orbits,'results':receipts,'scope':'Numerical diagnostic enumerating all1106 neutral nonvalence charge patterns, not implicit discovery or a rational certificate. For a strictly positive charge metric, worst Hubbard spin rows are obtained by alternating single occupations within every single run; global Sz=0 can be achieved because odd-run count is even. Relaxing metric positivity to nonnegative enlarges the candidate family. Negative optimized numerator floor would need exact dual reconstruction before an infeasibility claim.'}
Path('results/marginal_graded_hubbard8/charge_metric_diagnostic.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k!='feature_orbits'}))
