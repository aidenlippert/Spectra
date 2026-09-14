"""Cutting-plane discovery for finite-range positive charge products.

The charge-pattern oracle is the exact joint envelope DP.  This module never
materializes the 1106-pattern charge state space; it only retains witnesses
returned by the DP and transitions reachable from those witnesses.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse, json, math, time, sys, hashlib
import cvxpy as cp
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from experiments.marginal_joint_charge_dp import joint_envelope, replay

ROOT = Path(__file__).resolve().parents[1]

def factors(theta, keys, sites, radius):
    def f(x): return str(F(max(1, round(math.exp(float(x))*10**10)), 10**10))
    onsite = [[f(sum(theta[k]*v**a for k,(i,j,a,b) in enumerate(keys) if i==site and b==0)) for v in (-1,0,1)] for site in range(sites)]
    pairs=[]
    for i in range(sites):
        for j in range(i+1, min(sites, i+radius+1)):
            vals=[[f(sum(theta[k]*a**p*b**s for k,(ii,jj,p,s) in enumerate(keys) if ii==i and jj==j and s)) for b in (-1,0,1)] for a in (-1,0,1)]
            pairs.append({'i':i,'j':j,'values':vals})
    return onsite,pairs

def cert(onsite, pairs, data, gamma='-100', count_profile=None):
    return dict(data, kind='joint_charge_product_dp_v1',
                sites=data['sites'], modes=data['modes'], particles=data['particles'], U='4', t='1', target_lower=gamma,
                onsite=onsite, pairs=pairs, **({'doublon_weights':count_profile} if count_profile else {}))

def transitions(patterns, radius):
    edges=[]
    for r,q in enumerate(patterns):
        for i in range(len(q)-1):
            a,b=q[i],q[i+1]; targets=[]
            if abs(a-b)==1:
                z=list(q); z[i],z[i+1]=b,a; targets=[(z,1)]
            elif {a,b}=={-1,1}:
                z=list(q); z[i]=z[i+1]=0; targets=[(z,2)]
            elif a==b==0:
                for d in (-1,1):
                    z=list(q); z[i],z[i+1]=d,-d; targets.append((z,1))
            for z,w in targets:
                # Targets need not be retained as cuts, but their local
                # features are part of the source row's exact transition sum.
                edges.append((r,tuple(z),w))
    return edges

def run(radius, out, seconds=180, max_rounds=60, sites=8, target=None, resume_history=None, count_profile=False):
    if type(count_profile) is not bool or type(sites) is not int or sites%2 or not 2<=sites<=32 or radius not in (1,2,3):
        raise ValueError('Bounded even sites, local range and Boolean profile flag required')
    if not math.isfinite(seconds) or seconds<=0 or type(max_rounds) is not int or not 1<=max_rounds<=200:
        raise ValueError('Bounded positive search budgets required')
    started=time.monotonic(); out=Path(out)
    if out.exists() and any(out.iterdir()):raise ValueError('Preserve previous discovery artifacts')
    out.mkdir(parents=True,exist_ok=True)
    keys=[(i,i,a,0) for i in range(sites) for a in (1,2)]
    target_gamma=F(target) if target is not None else (F('-4.13') if radius==2 else F('-4.55'))
    profile=['0']+[str(i) for i in range(1,sites//2+1)]
    profile_counts=list(range(2,sites//2+1)) if count_profile else []
    keys += [(i,j,a,b) for i in range(sites) for j in range(i+1,min(sites,i+radius+1)) for a in (1,2) for b in (1,2)]
    # One initial cut per D, obtained from the DP oracle at uniform factors.
    from experiments.marginal_symbolic import mono, product as smul, scale, add, encode
    n=[mono(((1,i),(0,i))) for i in range(2*sites)]
    h=add(*(scale(smul(n[2*i],n[2*i+1]),4) for i in range(sites)),
          *(mono(((1,2*i+s),(0,2*(i+1)+s)),-1) for i in range(sites-1) for s in (0,1)),
          *(mono(((1,2*(i+1)+s),(0,2*i+s)),-1) for i in range(sites-1) for s in (0,1)))
    data={'sites':sites,'modes':2*sites,'particles':sites,'hamiltonian':encode(h)}
    one, pairs=factors(np.zeros(len(keys)),keys,sites,radius)
    seed=joint_envelope(cert(one,pairs,data,count_profile=profile),return_witnesses=True)
    patterns=[tuple(row['worst_charge_pattern']) for row in seed['rows']]
    source_history_sha=None
    if resume_history:
        raw=Path(resume_history).read_bytes()
        if len(raw)>50_000_000:raise ValueError('History byte budget exceeded')
        source_history_sha=hashlib.sha256(raw).hexdigest();old=json.loads(raw)
        if type(old) is not list or len(old)>200:raise ValueError('Invalid bounded history')
        for h in old:
            for row in h.get('dp_rows',[]):
                p=tuple(row.get('worst_charge_pattern',()))
                if len(p)!=sites or not all(type(x) is int and x in (-1,0,1) for x in p) or sum(p)!=0 or not any(p):raise ValueError('Invalid resumed physical charge cut')
                if p not in patterns:patterns.append(p)
    provenance={'sites':sites,'range':radius,'count_profile':count_profile,'target':str(target_gamma),
                'hamiltonian_sha256':hashlib.sha256(json.dumps({k:data[k] for k in ('modes','particles','hamiltonian')},sort_keys=True,separators=(',',':')).encode()).hexdigest(),
                'resume_history':str(resume_history) if resume_history else None,'resume_history_sha256':source_history_sha,
                'initial_patterns':len(patterns),'initial_metric':'Uniform local factors and g(D)=D',
                'generator_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                'scope':'Only bare uniform Hubbard input and DP-generated physical cuts; a resumed history supplies cuts only, never metric coefficients.'}
    (out/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    history=[]; accepted=None
    for iteration in range(max_rounds):
        if time.monotonic()-started > seconds: break
        edges=transitions(patterns,radius); q=np.asarray(patterns,dtype=int)
        D=np.sum(q*q,axis=1)/2
        X=np.array([q[:,i]**a*(q[:,j]**b if b else 1) for i,j,a,b in keys]+[(D==k).astype(float) for k in profile_counts]).T
        source=np.array([e[0] for e in edges],int)
        target_q=np.asarray([e[1] for e in edges],dtype=int)
        target_D=np.sum(target_q*target_q,axis=1)/2
        rates=np.array([e[2] for e in edges],float)*target_D/D[source]
        target_X=np.array([target_q[:,i]**a*(target_q[:,j]**b if b else 1) for i,j,a,b in keys]+[(target_D==k).astype(float) for k in profile_counts]).T
        nfeat=len(keys)+len(profile_counts)
        theta=cp.Variable(nfeat); eta=cp.Variable()
        rows=[]
        for r in range(len(patterns)):
            mask=source==r
            rows.append(cp.sum(cp.multiply(rates[mask],cp.exp((target_X[mask]-X[r])@theta))) if mask.any() else 0)
        prob=cp.Problem(cp.Minimize(eta),[cp.hstack(rows)-4*D<=eta,theta>=-10,theta<=10])
        remaining=max(0.1, seconds-(time.monotonic()-started))
        prob.solve(solver='CLARABEL',max_iter=150,time_limit=remaining,tol_gap_abs=1e-8,tol_feas=1e-8)
        if theta.value is None: raise RuntimeError('convex restricted problem failed: '+prob.status)
        onsite,pairs=factors(theta.value,keys,sites,radius)
        # Profile coefficients are appended after the local feature block.
        prof=['0','1']+[str(F(max(1,round(math.exp(float(theta.value[len(keys)+k-2]))*k*10**10)),10**10)) for k in profile_counts] if count_profile else profile
        (out/'proposal_checkpoint.json').write_text(json.dumps(dict(provenance,onsite=onsite,pairs=pairs,doublon_weights=prof,active_patterns=patterns),indent=2)+'\n')
        exact=joint_envelope(cert(onsite,pairs,data,count_profile=prof),return_witnesses=True)
        minimum=F(exact['minimum_lower']); added=0
        for row in exact['rows']:
            p=tuple(row['worst_charge_pattern'])
            if p not in patterns: patterns.append(p); added+=1
        history.append({'round':iteration,'restricted_patterns':len(patterns)-added,'generated_pattern_count':len(patterns),'added':added,'objective':float(prob.value),'dp_minimum':str(minimum),'dp_rows':exact['rows']})
        history_path=out/'history.json'; history_path.write_text(json.dumps(history,indent=2)+'\n')
        if minimum >= target_gamma:
            gamma=target_gamma
            final=dict(cert(onsite,pairs,data,str(gamma),count_profile=prof)); receipt=replay(final)
            (out/'certificate.json').write_text(json.dumps(final,indent=2)+'\n'); (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n'); accepted=receipt; break
        if not added: break
    (out/'history.json').write_text(json.dumps(history,indent=2)+'\n')
    (out/'numerical.json').write_text(json.dumps(dict(provenance,parameters=len(keys)+len(profile_counts),generated_pattern_count=len(patterns),rounds=len(history),accepted=accepted is not None,seconds=time.monotonic()-started,full_charge_enumeration=False),indent=2)+'\n')
    if accepted is None: raise RuntimeError('bounded cutting-plane run did not reach target')

if __name__=='__main__':
    p=argparse.ArgumentParser(); p.add_argument('--range',type=int,choices=(1,2,3),required=True); p.add_argument('--sites',type=int,default=8); p.add_argument('--target'); p.add_argument('--count-profile',action='store_true'); p.add_argument('--resume-history'); p.add_argument('--out',required=True); p.add_argument('--seconds',type=float,default=180); p.add_argument('--max-rounds',type=int,default=60); a=p.parse_args(); run(a.range,a.out,a.seconds,max_rounds=a.max_rounds,sites=a.sites,target=a.target,resume_history=a.resume_history,count_profile=a.count_profile)
