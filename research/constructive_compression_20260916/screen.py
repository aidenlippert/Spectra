"""Cold charge-MPS filtering screen. Numerical diagnostics are not certificates."""
import argparse
from fractions import Fraction as F
import json
from math import ceil,sqrt
from pathlib import Path
import time
import numpy as np
from research.constructive_compression_20260916 import charge_mps as cm,mps_numeric as mn
from research.constructive_compression_20260916.model import fixture,mpo,numerical_mpo,trace_seed


def solve(out,bond=32,rungs=4,steps=420,cycles=5,warm=0):
    if out.exists():raise FileExistsError(out)
    out.mkdir(parents=True);start=time.monotonic()
    op=mpo(rungs);W,wq=numerical_mpo(op);seed=trace_seed(rungs)
    b=float(8*rungs+2*(3*rungs-2))
    def h(state):return cm.apply_mpo(*state,W,wq)
    def energy(state):return mn.inner(state[0],h(state)[0])/mn.inner(state[0],state[0])
    def combine(states,coeff):return cm.linear_combination(states,coeff)
    def normed(state):
        a,q=state;a=[x.copy() for x in a];a[0]/=mn.norm(a);return a,q
    def step(prev,older,ell,first=False):
        alpha=(b+ell)/(b-ell);beta=-2/(b-ell);hp=h(prev)
        raw=combine([prev,hp] if first else [prev,hp,older],
                    [alpha,beta] if first else [2*alpha,2*beta,-1])
        a,q,d=cm.compress(*raw,bond)
        return (a,q),sqrt(d['discarded_squared_norm'])
    state=normed(seed);history=[]
    for cycle in range(cycles):
        e=energy(state);ell=e+.1
        old=state;prev,eta=step(state,None,ell,True)
        for j in range(2,21):
            new,eta=step(prev,old,ell);old,prev=prev,new
        state=normed(prev);e=energy(state)
        item={'cycle':cycle,'energy':e,'max_bond':max(len(q) for q in state[1]),
              'nonzero_entries':sum(np.count_nonzero(a) for a in state[0])}
        history.append(item);print(json.dumps(item),flush=True)
    ell=F(ceil(e*10**7)+1,10**7)
    (out/'upper_proposal.json').write_text(json.dumps({'charges':state[1],
                  'arrays':[a.tolist() for a in state[0]],'energy_float':e,'enumerated_configurations':0})+'\n')
    filter_seed=seed;overlap=1.;warm_history=[]
    for j in range(warm):
        hp=h(filter_seed)
        raw=combine([filter_seed,hp],[b/(b-float(ell)),-1/(b-float(ell))])
        a,q,dd=cm.compress(*raw,bond);eta=sqrt(dd['discarded_squared_norm']);nn=mn.norm(a)
        overlap=(overlap-eta)/nn;a[0]/=nn;filter_seed=(a,q)
        if (j+1)%50==0 or j+1==warm:
            ee=energy(filter_seed)
            item={'warm_step':j+1,'energy':ee,'uncertified_overlap_lower':overlap,'eta_float':eta}
            warm_history.append(item);print(json.dumps(item),flush=True)
        if overlap<=0:break
    if warm and overlap>0:
        e=energy(filter_seed);ell=F(ceil(e*10**7)+1,10**7)
        (out/'upper_proposal.json').write_text(json.dumps({'charges':filter_seed[1],
                    'arrays':[a.tolist() for a in filter_seed[0]],'energy_float':e,'enumerated_configurations':0})+'\n')
    delta=.0009;theta=np.arccosh(1+2*delta/(b-float(ell)))
    z=F(round(float(np.exp(-theta))*10**9),10**9);zz=float(z)
    lower=(F(str(b))+ell)/2-(F(str(b))-ell)*(z+1/z)/4
    old=filter_seed;prev,eta=step(filter_seed,None,float(ell),True);eta_sum=zz*eta
    stats=[]
    for j in range(1,steps+1):
        if j>1:
            new,eta=step(prev,old,float(ell));old,prev=prev,new;eta_sum+=zz**j*eta
        if j%20==0 or j==steps:
            vnorm=mn.norm(prev[0]);score=2*zz**j*vnorm+2*eta_sum/(1-zz*zz)
            item={'step':j,'norm_float':vnorm,'weighted_truncation_sum_float':eta_sum,
                  'uncertified_score':score,'uncertified_overlap_lower':overlap,'elapsed_seconds':time.monotonic()-start,
                  'max_bond':max(len(q) for q in prev[1]),'nonzero_entries':sum(np.count_nonzero(a) for a in prev[0])}
            stats.append(item);print(json.dumps(item),flush=True)
            if score<.9*overlap:break
            if 2*eta_sum/(1-zz*zz)>2*overlap:
                print('Stopping: numerical truncation contribution already exceeds criterion.',flush=True);break
    report={'status':'numerical_screen_only','bond_cap':bond,'rungs':rungs,
            'cold_history':history,'warm_history':warm_history,'chebyshev_history':stats,'ell':str(ell),'b':str(F(str(b))),
            'z':str(z),'proposed_lower':str(lower),'proposed_width_float':e-float(lower),
            'enumerated_configurations':0,'global_matrix_entries':0,
            'successful_numerical_score_only':stats[-1]['uncertified_score']<overlap,
            'elapsed_seconds':time.monotonic()-start}
    (out/'screen.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--bond',type=int,default=32)
    p.add_argument('--rungs',type=int,default=4);p.add_argument('--steps',type=int,default=420);p.add_argument('--cycles',type=int,default=5)
    p.add_argument('--warm',type=int,default=0)
    a=p.parse_args();solve(a.out,a.bond,a.rungs,a.steps,a.cycles,a.warm)
