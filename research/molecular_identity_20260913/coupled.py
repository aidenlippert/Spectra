"""Four-orbital negative directions coupled to the existing charge blocks."""
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import argparse
import hashlib
import json
import time

from experiments.marginal_symbolic import canonical,decode,encode,word_product,scale
from experiments.marginal_hunt_car import adj
from research.molecular_identity_20260913.local_blocks import local_groups

POLICY={'version':1,'seed':'creator_channels','candidate_width':4,
        'negative_directions_per_round':8,'enrichment_rounds':3,
        'charge_searched':-1,'adjoin_positive_charge':True,
        'coefficient_rounding_denominator':1000000,'minimum_violation':1e-7,
        'couple_to_existing_same_charge_symmetry_block':True,
        'solver':'CLARABEL','map_backend':'contraction','row_condition':False,
        'solver_seconds_per_solve':60,'local_particle_number_assumed':False}


def normalize(p):
    if not p:return None
    c=p[min(p)];return tuple(sorted(scale(p,1/c).items()))


def price(raw,h,m,signature,seen):
    import numpy as np
    started=time.monotonic()
    if raw['modes']!=m or decode(raw['hamiltonian'],m,4)!=h:
        raise ValueError('Pricing Hamiltonian differs')
    if len(raw['rows'])!=len(raw['dual']):raise ValueError('Moment length mismatch')
    y={tuple(tuple(letter) for letter in w):float(v) for w,v in zip(raw['rows'],raw['dual'])}
    if len(y)!=len(raw['rows']) or not all(np.isfinite(v) for v in y.values()):raise ValueError('Invalid numerical moments')
    candidates=[];pairs=0;blocks=0
    for support in combinations(range(m),4):
        for g in local_groups([support],signature):
            polys=g['polynomials']
            if sum(2*f-1 for f,_ in next(iter(polys[0])))!=-1:continue
            ws=[next(iter(p.items())) for p in polys];k=len(ws);matrix=np.zeros((k,k));blocks+=1
            for i,(w,c) in enumerate(ws):
                for j in range(i,k):
                    v,d=ws[j]
                    z=sum(float(t)*y.get(a,0.) for a,t in word_product(tuple((1-f,r) for f,r in reversed(w)),v))*float(c*d)
                    rz=sum(float(t)*y.get(a,0.) for a,t in word_product(tuple((1-f,r) for f,r in reversed(v)),w))*float(c*d)
                    matrix[i,j]=matrix[j,i]=(z+rz)/2;pairs+=2
            values,vectors=np.linalg.eigh(matrix)
            if values[0]>=-POLICY['minimum_violation']:continue
            coeffs=[F(round(float(v)*POLICY['coefficient_rounding_denominator']),POLICY['coefficient_rounding_denominator']) for v in vectors[:,0]]
            # Linear words are already present in the recipient block. Removing
            # them preserves the enlarged span and avoids redundant generators.
            p={w:c*z for (w,c),z in zip(ws,coeffs) if len(w)==3 and z}
            key=normalize(p)
            if key is None or key in seen:continue
            candidates.append({'support':support,'numeric_score':float(values[0]),'operator':p,'key':key})
    candidates.sort(key=lambda x:(x['numeric_score'],x['support'],x['key']))
    selected=[];new=set(seen)
    for item in candidates:
        if item['key'] in new:continue
        selected.append(item);new.add(item['key'])
        if len(selected)>=POLICY['negative_directions_per_round']:break
    return selected,{'supports_priced':len(list(combinations(range(m),4))),
                     'blocks_priced':blocks,'monomial_word_pair_products':pairs,
                     'seconds':time.monotonic()-started,
                     'scope':'Numerical pricing proposes exact rounded generators; it certifies no energy or omitted-cone bound.'}


def augment(base,signature,operators):
    parts={}
    for p in operators:
        for q in (p,canonical(adj(p))):
            w=next(iter(q));key=(sum(2*f-1 for f,_ in w),signature(w))
            parts.setdefault(key,[]).append(q)
    result=[]
    for (charge,sig),polys in sorted(parts.items()):
        target=next((g for g in base if sum(2*f-1 for f,_ in next(iter(g['polynomials'][0])))==charge
                    and signature(next(iter(g['polynomials'][0])))==sig),None)
        group={'name':f'coupled:{charge}:{sig}','polynomials':polys}
        if target is not None:group['merge_into']=target['name']
        result.append(group)
    return result


def run(fixture,out,resume=False):
    from research.certificate_scaling.commutator_dictionary import polynomial_groups,run as solve
    start=time.monotonic()
    if not resume:out.mkdir(parents=True,exist_ok=False)
    raw=fixture.read_bytes();f=json.loads(raw);m=f['modes'];h=decode(f['hamiltonian'],m,4)
    base,sig=polynomial_groups(h,m,creator_channels=True)
    seen={normalize(p) for g in base for p in g['polynomials']};operators=[];history=[]
    prior_elapsed=0.
    if resume:
        if json.loads((out/'policy.json').read_text())!=POLICY:raise ValueError('Resume policy differs')
        if (out/'experiment.json').exists():raise ValueError('Experiment already completed')
        history=json.loads((out/'history.json').read_text())
        failure=json.loads((out/'failure.json').read_text())
        prior_elapsed=failure['elapsed_seconds_from_file_timestamps']
    else:(out/'policy.json').write_text(json.dumps(POLICY,indent=2)+'\n')
    current=out/'seed'
    if resume:
        seed=json.loads((current/'certificate.json').read_text())
        if seed['modes']!=m or seed['particles']!=f['particles'] or decode(seed['hamiltonian'],m,4)!=h:
            raise ValueError('Resume Hamiltonian/sector differs')
    else:solve(h,m,f['particles'],current,creator_channels=True,map_backend='contraction',solver_seconds=POLICY['solver_seconds_per_solve'])
    for rnd in range(POLICY['enrichment_rounds']):
        if rnd<len(history):
            record=history[rnd]
            if record['round']!=rnd+1:raise ValueError('Resume round order differs')
            proposal_raw=(current/'dual_proposal.json').read_bytes()
            if hashlib.sha256(proposal_raw).hexdigest()!=record['source_proposal_sha256']:raise ValueError('Resume proposal changed')
            for item in record['selected']:
                p=decode(item['operator'],m,3);operators.append(p);seen.add(normalize(p))
        else:
            proposal_raw=(current/'dual_proposal.json').read_bytes();proposal=json.loads(proposal_raw)
            if proposal['particles']!=f['particles']:raise ValueError('Pricing particle number differs')
            selected,stats=price(proposal,h,m,sig,seen)
            if not selected:break
            for item in selected:operators.append(item['operator']);seen.add(item['key'])
            record={'round':rnd+1,'pricing':stats,
                    'selected':[{'support':x['support'],'numeric_score':x['numeric_score'],'operator':encode(x['operator'])} for x in selected],
                    'source_proposal_sha256':hashlib.sha256(proposal_raw).hexdigest()}
            history.append(record);(out/'history.json').write_text(json.dumps(history,indent=2)+'\n')
            print(json.dumps({'round':rnd+1,'added_generators_including_adjoints':2*len(selected),'pricing':stats}),flush=True)
        current=out/f'round_{rnd+1}'
        if (current/'receipt.json').exists():continue
        r=solve(h,m,f['particles'],current,creator_channels=True,map_backend='contraction',
                solver_seconds=POLICY['solver_seconds_per_solve'],additional_groups=augment(base,sig,operators))
        record.update({'receipt_path':str(current/'receipt.json'),'exact_lower':r['exact']['lower']})
        (out/'history.json').write_text(json.dumps(history,indent=2)+'\n')
    result={'policy':POLICY,'fixture_sha256':hashlib.sha256(raw).hexdigest(),
            'total_wall_seconds':time.monotonic()-start+prior_elapsed,'retained_generators_including_adjoints':2*len(operators),
            'resumed_after_channel_fix':resume,'prior_wall_seconds_estimated':prior_elapsed,
            'history':history,'many_body_space_enumerated':False,'source_factors_used':False,'source_upper_used':False,
            'scope':'Adaptive directions discovered from H and numerical duals; all costs included. Exact export determines accepted lower bounds.'}
    (out/'experiment.json').write_text(json.dumps(result,indent=2)+'\n')
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--fixture',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--resume',action='store_true')
    a=p.parse_args();run(a.fixture,a.out,a.resume)
