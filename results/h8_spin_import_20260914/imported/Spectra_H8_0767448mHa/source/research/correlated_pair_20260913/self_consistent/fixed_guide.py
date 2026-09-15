"""Fixed-guide adaptation of the paired self-consistent SOS construction.

The numerical proposal rounds coefficients, but acceptance reconstructs the
actual weighted squares. Constants are tracked in the bound, not iterated as
correction monomials. No expectation value replaces an operator.
"""
import argparse
from fractions import Fraction as F
import json
from pathlib import Path
import time
from experiments.marginal_symbolic import decode,encode,canonical,mono,add,scale,product,adj,hermitian,number_shift
from research.joint_patterns_20260913.core import anticommutator
from research.molecular_collective_20260913.core import digest


def ph(poly,n):
    return canonical({tuple((1-c if i<n else c,i) for c,i in w):v for w,v in poly.items()})


def guide(data):
    m,n=data['modes'],data['particles'];h=decode(data['hamiltonian'],m,4);hp=ph(h,n)
    diag=[hp.get(((1,i),(0,i)),F(0)) for i in range(m)]
    if not 0<n<m: raise ValueError('Both particle and hole modes required')
    lo=max(-diag[i] for i in range(n,m));hi=min(diag[:n])
    if lo>=hi: raise ValueError('No chemical potential gives a positive PH quadratic guide')
    mu=(lo+hi)/2;weights=[diag[i]+(-mu if i<n else mu) for i in range(m)]
    shifted=add(hp,scale(ph(number_shift(m,n),n),mu))
    base=add(*(mono(((1,i),(0,i)),w) for i,w in enumerate(weights)))
    V=add(shifted,scale(base,-1));V.pop((),None)
    return mu,weights,canonical(V),{'minimum_weight_Ha':str(min(weights)),'chemical_interval_Ha':[str(lo),str(hi)]}


def tau_map(V,weights):
    if not hermitian(V) or () in V or any(len(w) not in (2,4) for w in V): raise ValueError('Hermitian even nonscalar perturbation required')
    if any(F(v)<=0 for v in weights): raise ValueError('Positive guide weights required')
    out=[{} for _ in weights]
    for word,v in V.items():
        creators=[(p,i) for p,(c,i) in enumerate(word) if c];r=len(creators)
        if 2*r<len(word): continue # Its exact adjoint is supplied by the cross term.
        den=sum((weights[i] for p,i in creators),F(0));pref=F(1,2) if 2*r==len(word) else F(1)
        if den<=0: raise ValueError('Nonpositive perturbative denominator')
        for pos,i in creators:
            rest=word[:pos]+word[pos+1:]
            out[i]=add(out[i],mono(rest,v*pref*((-1)**pos)/den))
    return list(map(canonical,out))


def cross(taus,weights):
    return canonical(add(*(scale(add(product(mono(((1,i),)),t),product(adj(t),mono(((0,i),)))),weights[i]) for i,t in enumerate(taus))))


def expand(taus,weights):
    base=add(*(mono(((1,i),(0,i)),w) for i,w in enumerate(weights)))
    correction=cross(taus,weights)
    induced=add(*(scale(anticommutator(t,t),weights[i]) for i,t in enumerate(taus) if t))
    return canonical(add(base,correction,induced)),canonical(induced)


def check(data,cert):
    start=time.monotonic();m,n=data['modes'],data['particles']
    if set(cert)!={'kind','fixture_sha256','chemical_potential','weights','taus'} or cert['kind']!='fixed_guide_paired_sos_v1' or cert['fixture_sha256']!=digest(data): raise ValueError('Structured certificate binding/schema')
    if type(cert['chemical_potential']) is not str or any(type(w) is not str for w in cert['weights']): raise ValueError('Exact rational scalars required')
    weights=list(map(F,cert['weights']));mu=F(cert['chemical_potential'])
    if len(weights)!=m or len(cert['taus'])!=m or any(w<=0 for w in weights):raise ValueError('Positive guide weights required')
    taus=[decode(t,m,3) for t in cert['taus']]
    for i,t in enumerate(taus):
        charge=1 if i<n else -1
        if any(len(w) not in (1,3) for w in t) or any(sum((2*c-1)*(-1 if j<n else 1) for c,j in w)!=charge for w in t): raise ValueError('Odd physical-charge homogeneous correction required')
    h=decode(data['hamiltonian'],m,4)
    if not hermitian(h) or any(sum(2*c-1 for c,i in w) for w in h):raise ValueError('Invalid physical Hamiltonian')
    sos,induced=expand(taus,weights)
    remainder=add(ph(h,n),scale(ph(number_shift(m,n),n),mu),scale(sos,-1))
    if not hermitian(remainder) or any(len(w)>4 for w in remainder):raise ValueError('Exact paired closure failed')
    b=remainder.pop((),F(0));eps=sum(map(abs,remainder.values()),F(0));L=b-eps
    return {'status':'certified_lower','lower_Ha':str(L),'lower_float_Ha':float(L),'base_Ha':str(b),
            'residual_norm_bound_Ha':str(eps),'positive_factor_count':m+sum(bool(t) for t in taus),
            'tau_coefficients':sum(map(len,taus)),'residual_terms':len(remainder),
            'residual_max_degree':max(map(len,remainder),default=0),'no_fixed_N_enumeration':True,
            'no_full_cubic_Gram_solve':True,'replay_seconds':time.monotonic()-start}


def rounded(poly,den):
    return {w:F(round(c*den),den) for w,c in poly.items() if round(c*den)}


def build(data,iterations=5,damping=F(1,2),den=10**7):
    start=time.monotonic();mu,weights,target,diag=guide(data);current=target;history=[];best=None;bestrec=None
    def candidate(taus,label):
        nonlocal best,bestrec
        cert={'kind':'fixed_guide_paired_sos_v1','fixture_sha256':digest(data),'chemical_potential':str(mu),'weights':list(map(str,weights)),'taus':list(map(encode,taus))}
        rec=check(data,cert);history.append({'iteration':label,**rec})
        if bestrec is None or F(rec['lower_Ha'])>F(bestrec['lower_Ha']):best,bestrec=cert,rec
    candidate([{} for _ in weights],'guide_only')
    targetnorm=sum(map(abs,target.values()),F(0))
    for it in range(iterations):
        raw=tau_map(current,weights)
        # Exact cross identity is checked before rounding, on the actual PH perturbation.
        if cross(raw,weights)!=current:raise AssertionError('Creator distribution cross identity failed')
        taus=[rounded(t,den) for t in raw];candidate(taus,it)
        _,induced=expand(taus,weights);induced.pop((),None)
        current=rounded(add(scale(current,1-damping),scale(add(target,scale(induced,-1)),damping)),den)
        if sum(map(abs,current.values()),F(0))>4*targetnorm:
            history.append({'status':'coefficient_growth_stop','iteration':it,'threshold_times_initial_L1':4});break
    return best,{'best':bestrec,'history':history,'guide':diag,'wall_seconds':time.monotonic()-start,
                 'method':'Fixed PH guide adaptation; no rotated-overlap Hastings generalized solve; exact L1 residual.',
                 'discovery_uses_correlated_moments':False,'teacher_assisted':False}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('fixture');p.add_argument('output');p.add_argument('--iterations',type=int,default=5);a=p.parse_args()
    data=json.loads(Path(a.fixture).read_text());out=Path(a.output);out.mkdir(parents=True,exist_ok=True)
    cert,rec=build(data,a.iterations);(out/'certificate.json').write_text(json.dumps(cert,indent=2)+'\n');(out/'receipt.json').write_text(json.dumps(rec,indent=2)+'\n')
    print(json.dumps(rec,indent=2))
