"""Finite convex log-product discovery; exact joint DP is the separate gate."""
from pathlib import Path
from fractions import Fraction as F
import argparse,json,time,math
import numpy as np
from scipy.sparse import csr_matrix
import cvxpy as cp
from results.marginal_graded_hubbard8.discovery.charge_only_limit import build
from experiments.marginal_joint_charge_dp import replay

ROOT=Path(__file__).resolve().parents[1]


def run(radius,out,seconds=120,count_profile=False):
    started=time.monotonic();out=Path(out);out.mkdir(parents=True,exist_ok=True)
    patterns,A=build();q=np.asarray(patterns,dtype=int);sites=q.shape[1]
    D=np.sum(q*q,axis=1)/2
    keys=[(i,i,a,0) for i in range(sites) for a in (1,2)]
    keys += [(i,j,a,b) for i in range(sites) for j in range(i+1,min(sites,i+radius+1)) for a in (1,2) for b in (1,2)]
    X=np.array([q[:,i]**a*(q[:,j]**b if b else 1) for i,j,a,b in keys]).T
    if count_profile:X=np.c_[X,np.array([D==k for k in range(2,sites//2+1)],dtype=float).T]
    edges=[(r,c,value) for (r,c),value in A.items() if r!=c]
    source=np.array([r for r,c,v in edges]);target=np.array([c for r,c,v in edges])
    rates=np.array([v for r,c,v in edges],float)*D[target]/D[source]
    delta=csr_matrix(X[target]-X[source])
    collect=csr_matrix((np.ones(len(edges)),(source,np.arange(len(edges)))),shape=(len(q),len(edges)))
    theta=cp.Variable(X.shape[1]);eta=cp.Variable()
    penalties=collect@cp.multiply(rates,cp.exp(delta@theta))
    problem=cp.Problem(cp.Minimize(eta),[penalties-4*D<=eta,theta>=-10,theta<=10])
    problem.solve(solver='CLARABEL',max_iter=150,time_limit=seconds,tol_gap_abs=1e-8,tol_feas=1e-8)
    if theta.value is None:raise ValueError('No convex product proposal: '+str(problem.status))
    coeff=np.asarray(theta.value);rows=4*D-np.asarray(collect@(rates*np.exp(delta@coeff))).reshape(-1)
    def factor(value):return str(F(max(1,round(math.exp(float(value))*10**10)),10**10))
    onsite=[];pairs=[]
    for i in range(sites):
        onsite.append([factor(sum(coeff[k]*v**a for k,(ii,j,a,b) in enumerate(keys) if ii==i and not b)) for v in (-1,0,1)])
    for i in range(sites):
        for j in range(i+1,min(sites,i+radius+1)):
            values=[[factor(sum(coeff[k]*a**p*b**s for k,(ii,jj,p,s) in enumerate(keys) if ii==i and jj==j and s)) for b in (-1,0,1)] for a in (-1,0,1)]
            pairs.append(dict(i=i,j=j,values=values))
    data=json.loads((ROOT/'hamiltonian.json').read_text())
    gamma=F(math.floor(float(min(rows))*100)-1,100)
    certificate=dict(data,kind='joint_charge_product_dp_v1',sites=sites,U='4',t='1',target_lower=str(gamma),onsite=onsite,pairs=pairs)
    if count_profile:certificate['doublon_weights']=['0','1']+[factor(math.log(k)+coeff[len(keys)+k-2]) for k in range(2,sites//2+1)]
    (out/'certificate.json').write_text(json.dumps(certificate,indent=2)+'\n')
    numerical={'solver_status':problem.status,'objective':float(problem.value),'numerical_minimum_row':float(min(rows)),
               'range':radius,'parameters':len(coeff),'count_profile':count_profile,'charge_patterns_in_discovery':len(q),'transitions':len(edges),'discovery_seconds':time.monotonic()-started,
               'scope':'Finite enumerated convex proposal. This discovery still uses all1106chargepatterns; joint DP replay separately certifies rounded factors.'}
    (out/'numerical.json').write_text(json.dumps(numerical,indent=2)+'\n');print(json.dumps(numerical),flush=True)
    try:receipt=replay(certificate)
    except ValueError as error:
        (out/'rejection.json').write_text(json.dumps({'error':str(error)},indent=2)+'\n');raise
    (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--range',type=int,choices=(1,2,3),required=True);p.add_argument('--out',required=True);p.add_argument('--seconds',type=float,default=120);p.add_argument('--count-profile',action='store_true');a=p.parse_args()
    run(a.range,a.out,a.seconds,a.count_profile)
