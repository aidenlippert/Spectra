"""Selected-configuration proposal. Numerical discovery never accepts a bound."""
import argparse
from fractions import Fraction as F
from math import comb
import json
import os
from pathlib import Path
import resource
import time

import numpy as np
from scipy.linalg import eigh

from experiments.marginal_symbolic import add, encode, mono
from research.intervention_reduction_20260916.exact import decode_model, digest, integer_action


def beam_configurations(state, budget):
    arrays=[]
    for i,edges in enumerate(state['tensors']):
        a=np.zeros((len(state['bond_charges'][i]),2,len(state['bond_charges'][i+1])))
        for l,s,r,v in edges:
            a[l,s,r]=v/state['denominator']
        arrays.append(a)
    right=[None]*(len(arrays)+1);right[-1]=np.ones((1,1))
    for i in reversed(range(len(arrays))):
        right[i]=sum(arrays[i][:,s,:]@right[i+1]@arrays[i][:,s,:].T for s in (0,1))
    beam=[(0,np.ones(1),float(right[0][0,0]))]
    prefix_count=0
    for i,a in enumerate(arrays):
        candidates=[]
        for bits,v,_ in beam:
            for s in (0,1):
                w=v@a[:,s,:]
                weight=float(w@right[i+1]@w)
                prefix_count+=1
                if weight>0:
                    candidates.append((bits+(s<<i),w,weight))
        candidates.sort(key=lambda x:(-x[2],x[0]))
        beam=candidates[:budget]
    return [x[0] for x in beam], {'prefixes_evaluated':prefix_count,
        'selected_probability_float':sum(x[2] for x in beam)/float(right[0][0,0]),
        'probability_is_diagnostic_only':True}


def control_polynomials(data):
    spatial=data['modes']//2
    p=max(0,min(spatial-2,data['particles']//2-1));q=p+1
    density=add(*(mono(((1,2*p+s),(0,2*p+s)),1) for s in (0,1)),
                *(mono(((1,2*q+s),(0,2*q+s)),-1) for s in (0,1)))
    hop=add(*(mono(((1,2*p+s),(0,2*q+s)),1) for s in (0,1)),
            *(mono(((1,2*q+s),(0,2*p+s)),1) for s in (0,1)))
    return [density,hop], [p,q]


def propose(case, output, budget, ranks, amplitude=F(1,100), selection='energy'):
    start=time.monotonic()
    output.mkdir(parents=True,exist_ok=False)
    data=json.loads((case/'fixture.json').read_text())
    source=json.loads((case/'mps/state.json').read_text())
    if source['fixture_sha256']!=digest(data):
        raise ValueError('Source MPS belongs to another fixture')
    states,beam_info=beam_configurations(source,budget)
    if not states:
        raise ValueError('Empty selected space')
    h=decode_model(data)
    den,action=integer_action(h,states)
    index={s:i for i,s in enumerate(states)}
    matrix=np.zeros((len(states),len(states)))
    for j,col in enumerate(action):
        for state,value in col.items():
            if state in index:
                matrix[index[state],j]=value/den
    if np.max(np.abs(matrix-matrix.T))>1e-12:
        raise ValueError('Numerical Hermiticity diagnostic failed')
    controls,sites=control_polynomials(data)
    max_rank=min(max(ranks),len(states))
    if selection == 'energy':
        eigenvalues,vectors=eigh(matrix,subset_by_index=[0,max_rank-1],driver='evr')
        order=list(range(max_rank))
    else:
        eigenvalues,vectors=eigh(matrix,driver='evr')
        adjacency=np.zeros_like(matrix)
        control_matrices=[]
        for control in controls:
            dc,ac=integer_action(control,states)
            bc=np.zeros_like(matrix)
            for j,col in enumerate(ac):
                for state,value in col.items():
                    if state in index:
                        bc[index[state],j]=value/dc
            adjacency+=float(amplitude)*np.abs(vectors.T@bc@vectors)
            control_matrices.append(bc)
        np.fill_diagonal(adjacency,0)
        weights=np.zeros(len(states));weights[0]=1
        if selection == 'response':
            score=adjacency@weights
        elif selection == 'reachable':
            score=weights.copy();term=weights.copy()
            for k in range(1,17):
                term=(10/k)*(adjacency@term);score+=term
        elif selection == 'snapshot':
            ground=vectors[:,0]
            snapshots=[]
            for u in (-float(amplitude),float(amplitude)):
                for v in (-float(amplitude),float(amplitude)):
                    e,z=eigh(matrix+u*control_matrices[0]+v*control_matrices[1],driver='evr')
                    weights=z.T@ground
                    for t in (2.5,5.,7.5,10.):
                        state=z@(np.exp(-1j*(e-eigenvalues[0])*t)*weights)
                        snapshots.extend([state.real-ground*(ground@state.real),state.imag-ground*(ground@state.imag)])
            modes,sv,_=np.linalg.svd(np.column_stack(snapshots),full_matrices=False)
            # Reorthogonalize against the initial vector; use it exactly as column zero.
            basis=[ground]
            for column in modes.T:
                x=column.copy()
                for _ in range(2):
                    for b in basis:
                        x-=b*(b@x)
                norm=np.linalg.norm(x)
                if norm>1e-10:
                    basis.append(x/norm)
                if len(basis)>=max_rank:
                    break
            if len(basis)<max_rank:
                raise ValueError('Insufficient snapshot rank for requested candidates')
            vectors=np.column_stack(basis)
            order=list(range(max_rank))
        else:
            raise ValueError('Unknown selection rule')
        if selection != 'snapshot':
            order=[0]+sorted(range(1,len(states)),key=lambda i:(-score[i],i))
            vectors=vectors[:,order[:max_rank]]
    den_v=10**8
    outputs=[]
    for rank in sorted(set(min(r,len(states)) for r in ranks)):
        integer_vectors=np.rint(vectors[:,:rank]*den_v).astype(np.int64).tolist()
        proposal={'kind':'integer_control_subspace_v1','fixture_sha256':digest(data),
            'configurations':states,'vectors':integer_vectors,'denominator':den_v,
            'controls':[{'operator':encode(poly),'amplitude_Ha':str(amplitude)} for poly in controls]}
        name=f'proposal_r{rank}.json'
        (output/name).write_text(json.dumps(proposal,separators=(',',':'))+'\n')
        outputs.append(name)
    na,nb=source['spin_counts'];s=data['modes']//2
    report={'case':str(case),'fixture_sha256':digest(data),'inherited_MPS_sha256':digest(source),
        'inherited_state_discovery_cost_included':False,'input_integral_generation_included':False,
        'selected_budget':budget,'selected_configurations':len(states),
        'balanced_magnetic_sector_dimension':comb(s,na)*comb(s,nb),
        'selected_space_exhausts_magnetic_sector':len(states)==comb(s,na)*comb(s,nb),
        'reached_configurations':len(set().union(*(set(col) for col in action))),
        'selected_matrix_entries':len(states)**2,'lowest_eigenvalue_diagnostic_Ha':float(eigenvalues[0]),
        'selection_rule':selection,'eigenvectors_computed':len(eigenvalues),
        'selected_eigenvector_indices':None if selection=='snapshot' else order[:max_rank],
        'snapshot_eigensystems_computed':4 if selection=='snapshot' else 0,
        'control_spatial_orbitals':sites,'proposals':outputs,**beam_info,
        'discovery_seconds':time.monotonic()-start,'peak_RSS_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        'BLAS_threads':{k:os.environ.get(k) for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS')}}
    (output/'discovery.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('case',type=Path);p.add_argument('output',type=Path)
    p.add_argument('--budget',type=int,default=64);p.add_argument('--ranks',type=int,nargs='+',default=[4,8,16])
    p.add_argument('--selection',choices=['energy','response','reachable','snapshot'],default='energy')
    a=p.parse_args();propose(a.case,a.output,a.budget,a.ranks,selection=a.selection)
