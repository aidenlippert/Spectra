from __future__ import annotations
raise RuntimeError('Rejected exploratory prototype: use self_consistent.fixed_guide. The original creator map and closure diagnosis were incorrect.')
import json
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[3]; sys.path.insert(0,str(ROOT))
from fractions import Fraction as F
from experiments.marginal_symbolic import decode, verify, canonical, mono, add, scale, product, adj
from research.joint_patterns_20260913.core import anticommutator

def tau_from_perturbation(V, weights):
    """Distribute a balanced even perturbation over creator-removal maps.

    Returns odd charge-minus-one polynomials tau_i.  Each canonical monomial
    is assigned only to creator labels, with denominator sum(weights[creator]);
    equal creation/annihilation degree receives the prescribed half factor.
    """
    if any(F(v) <= 0 for v in weights.values()): raise ValueError('positive PH weights required')
    taus={i:{} for i in weights}
    for word, coeff in V.items():
        if not coeff: continue
        creators=[p for c,p in word if c]
        r=len(creators); annih=sum(1 for c,_ in word if not c)
        if r<annih: raise ValueError('monomial has fewer creators than annihilators')
        if not creators: raise ValueError('constant/annihilator-only perturbation unsupported')
        denom=sum((F(weights[i]) for i in creators),F(0)); pref=F(1,4) if r==annih else F(1,2)
        for pos,i in enumerate(creators):
            rest=word[:next(q for q,(c,p) in enumerate(word) if c and p==i)] + word[next(q for q,(c,p) in enumerate(word) if c and p==i)+1:]
            di=mono(((0,i),)); trial=product(adj(di),mono(rest)); lead=trial.get(word,F(0))
            if not lead: raise ValueError('creator-removal map has zero leading coefficient')
            taus[i]=add(taus[i],scale(mono(rest),F(coeff)*pref*F(weights[i],1)/denom/lead))
    return {i:canonical(p) for i,p in taus.items() if p}

def tau_cross(taus, weights):
    out={}
    for i,tau in taus.items():
        d=mono(((0,i),)); out=add(out,scale(add(product(adj(d),tau),product(adj(tau),d)),F(weights[i])))
    return canonical(out)

def self_consistent_update(Vtarget, weights, iterations=8, damping=F(1,2)):
    """Bounded fixed-point coefficient update; refuses unsupported charges."""
    if not 1 <= iterations <= 8 or not (F(0)<damping<=F(1)): raise ValueError('invalid iteration budget')
    current=canonical(Vtarget); history=[]
    for _ in range(iterations):
        taus=tau_from_perturbation(current,weights)
        induced={}
        for i,tau in taus.items():
            induced=add(induced,scale(add(product(adj(tau),tau),product(tau,adj(tau))),F(weights[i])))
        # Keep the full exact polynomial; no clipping or omitted terms.
        current=canonical(add(scale(current,F(1)-damping),scale(add(Vtarget,scale(induced,F(-1))),damping)))
        history.append({'terms':len(current),'induced_terms':len(induced)})
    return current,history

def guide_certificate(fixture: dict, correlated_moments=None, denominator=10**6):
    m,n=fixture['modes'],fixture['particles']; decode(fixture['hamiltonian'],m,4)
    if not isinstance(denominator,int) or denominator<=0: raise ValueError('positive denominator required')
    order=list(range(m))
    if correlated_moments is not None:
        if not isinstance(correlated_moments,dict): raise ValueError('moment map required')
        order.sort(key=lambda i:(-float(correlated_moments.get(str(i),0)),i))
    blocks=[{'name':f'guide_{i}','words':[[[0,i]]],'factor':[[1]]} for i in order]
    cert={'modes':m,'particles':n,'hamiltonian':fixture['hamiltonian'],'b':'0',
          'number_multiplier':[],'denominator':denominator,'blocks':blocks,
          'construction':'hf_particle_hole_guide_v1','correlated_moments_used':correlated_moments is not None}
    return cert,verify(cert)

def paired_cubic_certificate(fixture: dict, denominator=10**6):
    """Build paired linear/mixed-cubic odd factors with exact CAR cancellation."""
    m,n=fixture['modes'],fixture['particles']; decode(fixture['hamiltonian'],m,4)
    if m < 3: raise ValueError('at least three modes required')
    blocks=[]
    for i in range(m):
        j=(i+1)%m; k=(i+2)%m
        linear=((0,i),); cubic=next(iter(canonical(mono(((1,j),(0,k),(0,i))))))
        words=[[[c,p] for c,p in linear],[[c,p] for c,p in cubic]]
        blocks.append({'name':f'paired_{i}','words':words,'factor':[[1,1]]})
        aw=[[[1-c,p] for c,p in reversed(linear)],[[1-c,p] for c,p in reversed(cubic)]]
        blocks.append({'name':f'paired_{i}_adj','words':aw,'factor':[[1,1]]})
        anticommutator(mono(linear),mono(cubic))
    cert={'modes':m,'particles':n,'hamiltonian':fixture['hamiltonian'],'b':'0',
          'number_multiplier':[],'denominator':denominator,'blocks':blocks,
          'construction':'paired_linear_mixed_cubic_anticommutator_v1'}
    return cert,verify(cert)

def build_from_path(path:Path,out:Path):
    fixture=json.loads(path.read_text()); cert,rec=guide_certificate(fixture)
    out.mkdir(parents=True,exist_ok=True)
    (out/'certificate.json').write_text(json.dumps(cert,separators=(',',':'))+'\n')
    (out/'receipt.json').write_text(json.dumps(rec,indent=2)+'\n')
    return rec
