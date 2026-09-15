"""Split quartic Gram blocks with conserved-charge stabilizer representations."""
from fractions import Fraction as F
import json
from pathlib import Path

import numpy as np

from experiments.marginal_adaptive import assemble_and_solve
from experiments.marginal_coefficient import export, symmetric_upper
from experiments.marginal_collective import hopping_polynomial, verify_interval
from experiments.marginal_dual_exact import full_word_blocks, component_charge
from experiments.marginal_orbit_quartic import symmetry_group, word_action


def split_blocks(modes=10,degree=4,indices=None):
    original=full_word_blocks(modes,degree);group=symmetry_group(modes)
    if indices is not None:
        indices=set(indices)
        if any(type(i) is not int or i<0 or i>=len(original) for i in indices):
            raise ValueError('Invalid dictionary index')
    rng=np.random.default_rng(418);result=[];receipt=[]
    for index,block in enumerate(original):
        if indices is not None and index not in indices:continue
        words=block['words'];q=component_charge(words[0],modes);lookup={w:i for i,w in enumerate(words)}
        matrix=np.zeros((len(words),len(words)));order=0
        for mapping in group:
            # A mode in flavor i maps to flavor mapping[i] mod m.
            if any(q[mapping[i]%(modes//2)]!=q[i] for i in range(modes//2)):continue
            weight=int(rng.integers(-1000,1001));permutation=np.zeros_like(matrix)
            for j,w in enumerate(words):
                image,sign=word_action(w,mapping)
                permutation[lookup[image],j]=sign
            matrix+=weight*(permutation+permutation.T);order+=1
        eigenvalues,vectors=np.linalg.eigh(matrix);clusters=[]
        for j,value in enumerate(eigenvalues):
            if not clusters or abs(value-eigenvalues[clusters[-1][0]])>1e-7:clusters.append([j])
            else:clusters[-1].append(j)
        for j,cluster in enumerate(clusters):
            result.append({'name':f'{index}:{j}','words':words,'basis_transform':vectors[:,cluster]})
        receipt.append({'original_dimension':len(words),'group_order':order,'subspaces':list(map(len,clusters))})
    return result,receipt


def run():
    modes=10;t=F(1,5);h=hopping_polynomial(modes,t)
    blocks,decomposition=split_blocks(modes)
    out=Path(__file__).resolve().parents[1]/'results/marginal_stabilizer';out.mkdir(exist_ok=True)
    sizes=[b['basis_transform'].shape[1] for b in blocks]
    structural={'blocks':len(blocks),'largest':max(sizes),'psd_scalars':sum(n*(n+1)//2 for n in sizes),'decomposition':decomposition}
    (out/'structure.json').write_text(json.dumps(structural,indent=2)+'\n')
    print(json.dumps({k:structural[k] for k in ('blocks','largest','psd_scalars')}),flush=True)
    proposal,data=assemble_and_solve(h,modes,5,blocks,degree=4,groups=[[i,i+5] for i in range(5)],residual_penalty=True)
    cert,receipt=export(h,modes,5,blocks,proposal,denominator=10**9,operator_degree=4)
    cert['variational_upper']=symmetric_upper(modes,t);receipt.update(verify_interval(cert))
    receipt['numerical_objective']=proposal['numerical_objective']
    receipt['residual_penalty']=True
    (out/'certificate.json').write_text(json.dumps(cert)+'\n');(out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:receipt[k] for k in ('lower_float','width_float','status','build_seconds','solve_seconds')}),flush=True)
    return receipt


if __name__=='__main__':run()
