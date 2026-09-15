"""Seed a localized-basis signed-atom search with an exact empty DD proof."""
from fractions import Fraction as F
from itertools import combinations
import json
from pathlib import Path

from experiments.marginal_spin_reduction import SpinZeroOracle
from experiments.marginal_spin_constructor import spin_states, spin_blocks
from experiments.marginal_h6_complement import checked_blocks
from experiments.marginal_clique_gap import balanced_signs
from experiments.marginal_signed_atoms import verify_block, encode_vector, residual

def _weights(a):
    import numpy as np
    x=np.array([[float(v) for v in row] for row in a])
    c=np.diag(np.diag(x))-np.abs(x-np.diag(np.diag(x)))
    v=np.linalg.eigh(c)[1][:,0]
    w=np.maximum(1,np.rint(np.abs(v)*10**8).astype(int))
    return [int(z) for z in w]

def _dictionary(a):
    import numpy as np
    x=np.array([[abs(float(v)) for v in row] for row in a]); n=len(a); out=[]; seen=set()
    for i in range(n):
        top=sorted((j for j in range(n) if j!=i),key=lambda j:(-x[i,j],j))[:6]
        for k in (3,4):
            for tail in combinations(top,k-1):
                g=tuple(sorted((i,)+tail)); signs=balanced_signs(a,g)
                if signs is None: continue
                key=(g,tuple(signs))
                if key not in seen: seen.add(key); out.append(encode_vector(key))
    return out

def construct(source, output, retained_states, target=F(-6323,1000)):
    source, output=Path(source),Path(output)
    if output.exists(): raise ValueError('Preserve previous localized search seed')
    d=json.loads(source.read_text()); base={k:d[k] for k in ('modes','particles','hamiltonian')}
    oracle=SpinZeroOracle(base)
    if type(retained_states) is not list or len(retained_states)>32 or len(set(retained_states))!=len(retained_states): raise ValueError('Invalid retained states')
    oracle.retained(retained_states)
    groups=spin_blocks(oracle,retained_states); mats=checked_blocks(oracle,retained_states,[{'states':g} for g in groups],max_block_dimension=384)
    blocks=[]; dictionaries=[]; bounds=[]
    for a,g in zip(mats,groups):
        w=_weights(a); item={'states':g,'metric_weights':w,'atoms':[],'atom_scale':1}
        r,m=residual(a,w,[],1,True)
        baseline=min((r[i][i]-sum((abs(r[i][j]) for j in range(len(r)) if j!=i),F(0)))/m[i] for i in range(len(r)))
        verify_block(a,item,baseline,allow_amplitudes=True)
        blocks.append(item); bounds.append(baseline)
        dictionaries.append({'states':g,'metric_weights':w,'atom_dictionary':_dictionary(a)})
    cert=dict(base,retained_states=retained_states,kind='spin_rational_atom_complement_v1',target_lower=str(min(bounds)),blocks=blocks)
    # Verify the exact empty decomposition before writing any artifact.
    from experiments.marginal_signed_atoms import replay
    replay(cert)
    metric={'modes':d['modes'],'particles':d['particles'],'hamiltonian':d['hamiltonian'],'retained_states':retained_states,
            'kind':'spin_clique_metric_probe_v1','target_lower':str(F(target)),'blocks':blocks,
            'scope':'Localized exploratory metric and empty exact DD seed; finite H6 proposal.'}
    state=dict(metric,kind='spin_atom_dictionary_v1',vector_mode='rational',blocks=dictionaries,scope='Bounded initial direction dictionary; no global pricing claim.')
    output.mkdir(parents=True)
    for name,obj in [('metric_source.json',metric),('certificate.json',cert),('search_state.json',state)]: (output/name).write_text(json.dumps(obj,indent=2)+'\n')
    return {'target_lower':str(min(bounds)),'requested_target_lower':str(F(target)),
            'requested_target_certified':min(bounds)>=F(target),
            'blocks':len(blocks),'dictionary_atoms':sum(len(x['atom_dictionary']) for x in dictionaries),'retained_states':retained_states}

if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--source',required=True);p.add_argument('--output',required=True);p.add_argument('--selection-source',required=True);p.add_argument('--target',default='-6.323')
    a=p.parse_args(); retained=json.loads(Path(a.selection_source).read_text())['localized']['p_states']
    print(json.dumps(construct(a.source,a.output,retained,F(a.target)),indent=2))
