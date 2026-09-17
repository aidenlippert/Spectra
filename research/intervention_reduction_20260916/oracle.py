"""Numerical full magnetic-sector propagation for independent diagnostics.

This explicitly enumerates its reference sector and is never imported by an
accepting checker. It consumes raw supplied CAR words with a separate ladder
implementation and uses sparse exponential propagation.
"""
import argparse
from fractions import Fraction as F
from itertools import combinations
import json
from pathlib import Path
import resource
import time
import numpy as np
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import expm_multiply
from research.intervention_reduction_20260916.exact import digest


def literal_action(word,label,modes):
    occupied=[p for p in range(modes) if (label>>p)&1]
    sign=1
    for creation,p in reversed(word):
        found=p in occupied
        if (creation and found) or (not creation and not found):
            return None
        below=sum(q<p for q in occupied)
        if below%2:
            sign=-sign
        if creation:
            occupied.append(p)
        else:
            occupied.remove(p)
    return sum(1<<p for p in occupied),sign


def matrix(terms,labels,modes):
    index={s:i for i,s in enumerate(labels)};entries=[];rows=[];ptr=[0]
    for s in labels:
        column={}
        for term in terms:
            word=term['word']
            if any(sum((2*c-1) for c,p in word if p%2==spin) for spin in (0,1)):
                raise ValueError('Numerical oracle requires separate spin-count conservation')
            result=literal_action(word,s,modes)
            if result is not None:
                target,sign=result
                column[index[target]]=column.get(index[target],0.)+sign*float(F(term['coefficient']))
        for row,value in column.items():
            if value:
                rows.append(row);entries.append(value)
        ptr.append(len(entries))
    return csc_matrix((entries,rows,ptr),shape=(len(labels),len(labels)))


def run(data,proposal,trajectory,receipt):
    start=time.monotonic();m=data['modes'];s=m//2;first=proposal['configurations'][0]
    counts=[sum((first>>p)&1 for p in range(spin,m,2)) for spin in (0,1)]
    if any([sum((label>>p)&1 for p in range(spin,m,2)) for spin in (0,1)]!=counts for label in proposal['configurations']):
        raise ValueError('Oracle initial configurations have different spin counts')
    labels=[sum(1<<(2*p) for p in a)+sum(1<<(2*p+1) for p in b)
            for a in combinations(range(s),counts[0]) for b in combinations(range(s),counts[1])]
    idx={label:i for i,label in enumerate(labels)}
    v=np.zeros((len(labels),len(proposal['vectors'][0])))
    for label,row in zip(proposal['configurations'],proposal['vectors']):
        v[idx[label]]=np.asarray(row)/proposal['denominator']
    initial=v[:,0]/np.linalg.norm(v[:,0]);state=initial.astype(complex)
    operators=[matrix(data['hamiltonian'],labels,m)]+[matrix(c['operator'],labels,m) for c in proposal['controls']]
    built=time.monotonic();groups=[]
    for segment in trajectory['segments']:
        controls=tuple(segment['controls_Ha']);duration=F(segment['duration_atomic_time'])
        if groups and groups[-1][0]==controls:
            groups[-1]=(controls,groups[-1][1]+duration)
        else:
            groups.append((controls,duration))
    for controls,duration in groups:
        h=operators[0]+sum(float(F(u))*a for u,a in zip(controls,operators[1:]))
        state=expm_multiply((-1j*float(duration))*h,state,traceA=(-1j*float(duration))*h.diagonal().sum())
    endpoint=np.asarray(trajectory['segments'][-1]['coefficients'],dtype=float)
    c=(endpoint[:,:,0]+1j*endpoint[:,:,1]).sum(axis=0)/trajectory['coefficient_denominator']
    reduced=v@c/np.linalg.norm(v[:,0])
    # The proposal removes a fixed energy phase.
    aligned=state*np.exp(1j*float(F(trajectory['phase_shift_Ha']))*float(F(trajectory['horizon_atomic_time'])))
    delta=float(np.linalg.norm(aligned-reduced/np.linalg.norm(reduced)))
    orbital=trajectory['observed_spatial_orbital']
    occupations=np.array([((label>>(2*orbital))&1)+((label>>(2*orbital+1))&1) for label in labels])
    population=float(np.vdot(state,occupations*state).real)
    lo,hi=map(lambda x:float(F(x)),receipt['population_interval'])
    error_bound=float(F(receipt['normalized_state_error_bound']))
    return {'scope':'uncertified numerical diagnostic on explicitly enumerated full magnetic sector',
        'fixture_sha256':digest(data),'proposal_sha256':digest(proposal),'reference_dimension':len(labels),
        'reference_sparse_entries':[a.nnz for a in operators], 'construction_seconds':built-start,
        'propagation_and_comparison_seconds':time.monotonic()-built,'total_seconds':time.monotonic()-start,
        'reference_population':population,'certified_population_interval':[lo,hi],
        'observed_state_error':delta,'certified_state_error':error_bound,
        'inside_state_bound':delta<=error_bound+1e-10,'inside_population_interval':lo-1e-10<=population<=hi+1e-10,
        'reference_norm':float(np.linalg.norm(state)),'peak_RSS_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('fixture',type=Path);p.add_argument('proposal',type=Path)
    p.add_argument('trajectory',type=Path);p.add_argument('receipt',type=Path);p.add_argument('output',type=Path);a=p.parse_args()
    if a.output.exists():
        raise FileExistsError(a.output)
    result=run(*(json.loads(path.read_text()) for path in (a.fixture,a.proposal,a.trajectory,a.receipt)))
    a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)
    if not result['inside_state_bound'] or not result['inside_population_interval']:
        raise SystemExit('Numerical oracle contradicts accepted bound')
