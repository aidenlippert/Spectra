"""Price bounded four-orbital blocks without constructing a full cubic map."""
from itertools import combinations
from pathlib import Path
import argparse
import hashlib
import json
import time

from experiments.marginal_symbolic import decode,word_product
from research.molecular_identity_20260913.local_blocks import local_groups

POLICY={'version':1,'seed':'creator_channels_plus_all_three_orbital_blocks',
        'candidate_width':4,'new_supports_per_round':8,'enrichment_rounds':2,
        'score':'minimum_eigenvalue_across_local_charge_and_exact_symmetry_blocks',
        'solver':'CLARABEL','map_backend':'contraction','row_condition':False,
        'solver_seconds_per_solve':60,'local_particle_number_assumed':False}


def price(raw,h,m,signature,excluded=()):
    import numpy as np
    if raw['modes']!=m or decode(raw['hamiltonian'],m,4)!=h:
        raise ValueError('Pricing Hamiltonian differs')
    if len(raw['rows'])!=len(raw['dual']):raise ValueError('Moment length mismatch')
    y={tuple(tuple(letter) for letter in w):float(v) for w,v in zip(raw['rows'],raw['dual'])}
    if len(y)!=len(raw['rows']) or not all(np.isfinite(v) for v in y.values()):
        raise ValueError('Invalid numerical moments')
    proposals=[];pairs=0;blocks=0;started=time.monotonic()
    for support in combinations(range(m),4):
        if support in excluded:continue
        gs=local_groups([support],signature);minimum=0.
        for g in gs:
            # These candidates are monomials; no global coefficient map is made.
            ws=[next(iter(p.items())) for p in g['polynomials']];k=len(ws)
            matrix=np.zeros((k,k));blocks+=1
            for i,(w,c) in enumerate(ws):
                for j in range(i,k):
                    v,d=ws[j]
                    forward=sum(float(z)*y.get(a,0.) for a,z in word_product(tuple((1-f,t) for f,t in reversed(w)),v))*float(c*d)
                    reverse=sum(float(z)*y.get(a,0.) for a,z in word_product(tuple((1-f,t) for f,t in reversed(v)),w))*float(c*d)
                    matrix[i,j]=matrix[j,i]=(forward+reverse)/2;pairs+=2
            minimum=min(minimum,float(np.linalg.eigvalsh(matrix)[0]))
        proposals.append({'support':support,'numeric_score':minimum})
    proposals.sort(key=lambda x:(x['numeric_score'],x['support']))
    return proposals,{'supports_priced':len(proposals),'blocks_priced':blocks,
                      'monomial_word_pair_products':pairs,'seconds':time.monotonic()-started,
                      'scope':'Heuristic pricing only. Missing numerical moments are zero; no omitted-cone bound is asserted.'}


def run(fixture,out,seed=None):
    from research.certificate_scaling.adaptive_block_discovery import partition
    from research.certificate_scaling.commutator_dictionary import run as solve
    start=time.monotonic();out.mkdir(parents=True,exist_ok=False)
    raw=fixture.read_bytes();f=json.loads(raw);m=f['modes'];h=decode(f['hamiltonian'],m,4)
    _,sig,_=partition(h,m,'quadratic',True)
    (out/'policy.json').write_text(json.dumps(POLICY,indent=2)+'\n')
    supports=list(combinations(range(m),3));chosen=[];history=[]
    if seed is None:
        seed=out/'seed'
        solve(h,m,f['particles'],seed,solver_seconds=POLICY['solver_seconds_per_solve'],
              creator_channels=True,map_backend='contraction',additional_groups=local_groups(supports,sig))
    seed_receipt=json.loads((seed/'receipt.json').read_text())
    current=seed
    for rnd in range(POLICY['enrichment_rounds']):
        proposal_raw=(current/'dual_proposal.json').read_bytes();proposal=json.loads(proposal_raw)
        if proposal['particles']!=f['particles']:raise ValueError('Pricing particle number differs')
        ranked,stats=price(proposal,h,m,sig,set(chosen))
        selected=[tuple(row['support']) for row in ranked if row['numeric_score'] < -1e-7][:POLICY['new_supports_per_round']]
        if not selected:break
        chosen.extend(selected);groups=local_groups(supports+chosen,sig)
        current=out/f'round_{rnd+1}'
        record={'round':rnd+1,'selected_supports':selected,'pricing':stats,
                'top_candidates':ranked[:POLICY['new_supports_per_round']],
                'source_proposal_sha256':hashlib.sha256(proposal_raw).hexdigest()}
        history.append(record)
        (out/'history.json').write_text(json.dumps(history,indent=2)+'\n')
        print(json.dumps({'stage':'selected','round':rnd+1,'pricing':stats,'selected_supports':selected}),flush=True)
        r=solve(h,m,f['particles'],current,solver_seconds=POLICY['solver_seconds_per_solve'],
                creator_channels=True,map_backend='contraction',additional_groups=groups)
        record['receipt_path']=str(current/'receipt.json');record['exact_lower']=r['exact']['lower']
        (out/'history.json').write_text(json.dumps(history,indent=2)+'\n')
    result={'policy':POLICY,'fixture_sha256':hashlib.sha256(raw).hexdigest(),
            'seed':str(seed),'seed_receipt_sha256':hashlib.sha256((seed/'receipt.json').read_bytes()).hexdigest(),
            'seed_discovery_seconds':seed_receipt['wall_seconds'],
            'additional_wall_seconds':time.monotonic()-start,
            'seed_reused':seed.parent!=out,'selected_four_orbital_supports':chosen,
            'history':history,'many_body_space_enumerated':False,
            'source_factors_used':False,'source_upper_used':False,
            'scope':'Hamiltonian-only adaptive rule. Accuracy and runtime are empirical; exact exported lower bounds are checked separately.'}
    (out/'experiment.json').write_text(json.dumps(result,indent=2)+'\n')
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--fixture',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True);p.add_argument('--seed',type=Path)
    a=p.parse_args();run(a.fixture,a.out,a.seed)
