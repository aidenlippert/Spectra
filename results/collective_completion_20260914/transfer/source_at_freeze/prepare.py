"""Construct paired maps and guide moments directly at quartic degree.

No sextic coefficient table, three-particle moment matrix, or prior Gram
certificate is constructed. Full cubic word-pair maps remain polynomial
preparation work and are explicitly counted.
"""
import argparse
from fractions import Fraction as F
from itertools import combinations
import json
from pathlib import Path
import time
import numpy as np
from scipy import sparse
from experiments.marginal_symbolic import decode
from experiments.marginal_coefficient import gram_map,dagger
from research.certificate_scaling.adaptive_block_discovery import partition
from research.reconstruction_compression_20260914.inputs import dump,sha
from research.reconstruction_compression_20260914.moments import Proposal

def build(case,fixture,state,out):
    t=time.monotonic();folder=Path(out);folder.mkdir(parents=True,exist_ok=False)
    data=json.loads(Path(fixture).read_text());h=decode(data['hamiltonian'],data['modes'],4);m,n=data['modes'],data['particles']
    groups,sig,stats=partition(h,m,'mixed',True)
    if n>6:groups=[{**g,'words':[w for w in g['words'] if len(w)!=1]} if g['name'].startswith('mixed') else g for g in groups]
    groups=[g for g in groups if g['words']]
    rows=[tuple((1,i) for i in a)+tuple((0,i) for i in b) for k in range(3)
          for a in combinations(range(m),k) for b in combinations(range(m),k)
          if a<=b and sig(tuple((1,i) for i in a)+tuple((0,i) for i in b))==sig(())]
    lookup={w:i for i,w in enumerate(rows)};used=set();pairs=[];work=0
    for i,g in enumerate(groups):
        if i in used:continue
        words=g['words'];size=len(words);M=gram_map(words,lookup).tocsc();work+=size**2
        if max(map(len,words))>2:
            found=None
            for j,other in enumerate(groups):
                if j==i or j in used or len(other['words'])!=size:continue
                loc={w:k for k,w in enumerate(other['words'])}
                if all(dagger(w) in loc for w in words):found=j;order=np.array([loc[dagger(w)] for w in words]);break
            if found is None:raise ValueError('Odd dictionary lacks adjoint partner')
            j=found;other=gram_map(groups[j]['words'],lookup).tocsc();work+=size**2
            M+=other[:,(order[None,:]*size+order[:,None]).ravel()];used.add(j)
            pairs.append({'members':[i,j],'adjoint_order':order.tolist()})
        else:pairs.append({'members':[i]})
        M.eliminate_zeros();sparse.save_npz(folder/f'paired_map_{len(pairs)-1}.npz',M)
    map_end=time.monotonic();oracle=Proposal(data,json.loads(Path(state).read_text()));oracle.fill(rows)
    weight=np.array([1 if tuple(i for c,i in w if c)==tuple(i for c,i in w if not c) else 2 for w in rows])
    y=weight*np.array([oracle.cache[w] for w in rows]);Cs={}
    for k,pair in enumerate(pairs):
        size=len(groups[pair['members'][0]]['words']);M=sparse.load_npz(folder/f'paired_map_{k}.npz')
        C=np.asarray(M.T@y).reshape(size,size);Cs[f'C_{k}']=(C+C.T)/2
    np.savez_compressed(folder/'moments.npz',**Cs)
    meta={'kind':'paired_quartic_v1','case':case,'fixture':str(fixture),'state':str(state),
          'fixture_sha256':sha(fixture),'state_sha256':sha(state),'groups':groups,'rows':rows,'pairs':pairs,'stats':stats}
    dump(folder/'frame.json',meta)
    rec={'case':case,'map_seconds':map_end-t,'moment_seconds':time.monotonic()-map_end,'total_seconds':time.monotonic()-t,
         'coefficient_rows':len(rows),'maximum_moment_degree':max(map(len,rows)),'requested_moments':len(rows),
         'cached_moments':len(oracle.cache),'transfer_steps':oracle.steps,'cubic_word_pairs_visited':work,
         'map_nonzeros':sum(sparse.load_npz(folder/f'paired_map_{k}.npz').nnz for k in range(len(pairs))),
         'minimum_guide_eigenvalue':min(float(np.linalg.eigvalsh(C)[0]) for C in Cs.values()),
         'many_body_states_enumerated':0,'full_cubic_teacher_used':False,'three_particle_moments_constructed':False}
    dump(folder/'receipt.json',rec);print(json.dumps(rec),flush=True)

def select(meta,Cs,rank):
    groups=meta['groups'];maps=[];guide=[]
    for k,pair in enumerate(meta['pairs']):
        i=pair['members'][0];size=len(groups[i]['words'])
        if len(pair['members'])==1:V=np.eye(size);members=[(i,V)]
        else:
            ev,V=np.linalg.eigh(Cs[f'C_{k}']);V=np.rint(V[:,:min(size,rank)]*1e8)/1e8
            j=pair['members'][1];W=np.empty_like(V);W[pair['adjoint_order']]=V;members=[(i,V),(j,W)]
            guide.append({'name':groups[i]['name'],'lowest':float(ev[0]),'highest_retained':float(ev[min(size,rank)-1])})
        maps.append((members,groups[i]['name']))
    return maps,guide

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('case');p.add_argument('fixture');p.add_argument('state');p.add_argument('out');a=p.parse_args()
    build(a.case,a.fixture,a.state,a.out)
