"""Bounded numerical pulse design in a reused molecular subspace.

The exact robust replay, performed separately, decides success. Searching the
finite menu does not establish globally optimal control or laboratory reachability.
"""
import argparse
from fractions import Fraction as F
from itertools import product
import json
from pathlib import Path
import resource
import time
import numpy as np
from scipy.linalg import expm
from research.intervention_reduction_20260916.predict import projected_generators,propose
from research.intervention_reduction_20260916.exact import digest


def design(data,proposal,output):
    started=time.monotonic();output.mkdir(exist_ok=False,parents=True)
    v,aa=projected_generators(data,proposal)
    rank=v.shape[1];initial=np.zeros(rank,dtype=complex);initial[0]=1
    orbital=min(data['modes']//2-1,data['particles']//2)
    occ=np.array([((s>>(2*orbital))&1)+((s>>(2*orbital+1))&1) for s in proposal['configurations']])
    observable=v.T@(occ[:,None]*v);metric=v.T@v
    initial_population=observable[0,0]/metric[0,0]
    amplitudes=[F(c['amplitude_Ha']) for c in proposal['controls']]
    choices=list(product(*[[-a,F(0),a] for a in amplitudes]))
    durations=[3,4,3];shift=aa[0][0,0]
    propagators={(t,i):expm(-1j*t*(aa[0]-shift*np.eye(rank)+sum(float(u)*a for u,a in zip(control,aa[1:]))))
                 for t in set(durations) for i,control in enumerate(choices)}
    candidates=[]
    for indices in product(range(len(choices)),repeat=len(durations)):
        state=initial.copy()
        for t,i in zip(durations,indices):
            state=propagators[t,i]@state
        value=float((np.vdot(state,observable@state)/np.vdot(state,metric@state)).real)
        candidates.append((value-initial_population,indices))
    candidates.sort(key=lambda x:(-x[0],x[1]))
    exported=[]
    for k,(score,indices) in enumerate(candidates[:3],1):
        schedule=[(F(t),list(choices[i])) for t,i in zip(durations,indices)]
        trajectory,info=propose(data,proposal,schedule,order=18,stable=True)
        path=output/f'candidate_{k}.json'
        path.write_text(json.dumps(trajectory,separators=(',',':'))+'\n')
        exported.append({'file':path.name,'predicted_population_increase':score,**info})
    report={'scope':'finite menu search on existing reduced model; numerical scores only',
        'fixture_sha256':digest(data),'proposal_sha256':digest(proposal),'reduced_dimension':rank,
        'requested_population_increase':'1/40','requested_uniform_control_deviation_Ha':'1/400000',
        'amplitude_grid':'negative maximum, zero, positive maximum per control',
        'durations_atomic_time':durations,'candidates_evaluated':len(candidates),
        'discovery_repeated':False,'parent_discovery_cost_included':False,
        'top_candidates':exported,'total_seconds':time.monotonic()-started,
        'peak_RSS_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (output/'design.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('fixture',type=Path);p.add_argument('proposal',type=Path)
    p.add_argument('output',type=Path);a=p.parse_args()
    design(json.loads(a.fixture.read_text()),json.loads(a.proposal.read_text()),a.output)
