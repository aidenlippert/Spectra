"""Add collective directions from the actual driven residual, without a dense
Hamiltonian diagonalization on the enlarged configuration set.

This is a matrix-free residual/POD proposal, not an enumeration-free claim:
every reached configuration and all inherited discovery are charged.
"""
import argparse
from fractions import Fraction as F
import json
from pathlib import Path
import resource
import time
import numpy as np
from scipy.linalg import expm
from scipy.sparse import csc_matrix
from experiments.marginal_symbolic import decode
from research.intervention_reduction_20260916.exact import decode_model, digest, integer_action


def enrich(data,old,added):
    started=time.monotonic()
    states=old['configurations'];v=np.asarray(old['vectors'],dtype=float)/old['denominator']
    rank=v.shape[1];polys=[decode_model(data)]+[decode(c['operator'],data['modes'],4) for c in old['controls']]
    raw=[integer_action(p,states) for p in polys]
    labels=sorted(set(states).union(*(set(c) for _,cols in raw for c in cols)))
    index={s:i for i,s in enumerate(labels)}
    expanded=np.zeros((len(labels),rank))
    for s,row in zip(states,v):
        expanded[index[s]]=row
    acted=[]
    for den,columns in raw:
        entries=[];row_indices=[];pointers=[0]
        for column in columns:
            for s,value in column.items():
                entries.append(value/den);row_indices.append(index[s])
            pointers.append(len(entries))
        operator=csc_matrix((entries,row_indices,pointers),shape=(len(labels),len(states)))
        acted.append(operator@v)
    g=expanded.T@expanded
    reduced=[np.linalg.solve(g,expanded.T@a) for a in acted]
    residuals=[a-expanded@k for a,k in zip(acted,reduced)]
    radii=[float(F(c['amplitude_Ha'])) for c in old['controls']]
    if len(radii)!=2:
        raise ValueError('Frozen discovery rule uses two controls')
    e0=np.zeros(rank);e0[0]=1;shift=reduced[0][0,0]
    snapshots=[]
    for u in (-radii[0],radii[0]):
        for w in (-radii[1],radii[1]):
            a=reduced[0]+u*reduced[1]+w*reduced[2]-shift*np.eye(rank)
            r=residuals[0]+u*residuals[1]+w*residuals[2]
            for t in (0.,2.5,5.,7.5,10.):
                state=expm(-1j*t*a)@e0
                defect=r@state
                snapshots.extend([defect.real,defect.imag])
    modes,sv,_=np.linalg.svd(np.column_stack(snapshots),full_matrices=False)
    orthogonal,_=np.linalg.qr(expanded,mode='reduced')
    new=[]
    for vector in modes.T:
        x=vector.copy()
        for _ in range(2):
            x-=orthogonal@(orthogonal.T@x)
            for y in new:
                x-=y*(y@x)
        norm=np.linalg.norm(x)
        if norm>1e-10:
            new.append(x/norm)
        if len(new)==added:
            break
    if len(new)<added:
        raise ValueError('Insufficient independent residual directions')
    original_rows=dict(zip(states,old['vectors']))
    rounded=np.rint(np.column_stack(new)*old['denominator']).astype(np.int64).tolist()
    vectors=[original_rows.get(s,[0]*rank)+row for s,row in zip(labels,rounded)]
    proposal={**old,'configurations':labels,'vectors':vectors}
    return proposal,{'parent_sha256':digest(old),'proposal_sha256':digest(proposal),
        'old_rank':rank,'new_rank':rank+added,'old_configurations':len(states),'reached_configurations':len(labels),
        'configuration_matrix_diagonalized':False,'all_reached_labels_materialized':True,
        'training_controls':'unchanged four corners and times 0,2.5,5,7.5,10 au',
        'singular_values_diagnostic':sv.tolist(),'construction_seconds':time.monotonic()-started,
        'peak_RSS_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        'inherited_discovery_included':False}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('fixture',type=Path);p.add_argument('proposal',type=Path)
    p.add_argument('output',type=Path);p.add_argument('--added',type=int,default=8);a=p.parse_args()
    if a.output.exists():
        raise FileExistsError(a.output)
    result,report=enrich(json.loads(a.fixture.read_text()),json.loads(a.proposal.read_text()),a.added)
    a.output.write_text(json.dumps(result,separators=(',',':'))+'\n')
    a.output.with_suffix('.construction.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report),flush=True)
