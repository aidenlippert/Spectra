"""Numerically propose rational polynomial trajectories; exact replay is separate."""
import argparse
from fractions import Fraction as F
import json
from pathlib import Path
import time
import numpy as np
from research.intervention_reduction_20260916.exact import digest, integer_action, decode_model
from experiments.marginal_symbolic import decode


POLICIES={
    'heldout_switch':[(F(5),[F(1,100),F(0)]),(F(5),[F(0),F(-1,100)])],
    'heldout_interior':[(F(10),[F(3,1000),F(-7,1000)])],
    'heldout_three':[(F(3),[F(-1,200),F(1,125)]),(F(4),[F(1,100),F(-1,250)]),(F(3),[F(-1,250),F(-1,100)])],
}


def stable_step(generator, requested):
    """Numerical proposal safeguard; the residual checker certifies its result."""
    norm=float(np.linalg.norm(generator,ord=np.inf))
    # h ||A||_infinity <= 2. For order 18 the Taylor remainder is < 4e-11
    # per step. This is only a planning estimate; no floating norm is trusted.
    divisor=max(1,int(np.ceil(float(requested)*norm/2)))
    return requested/divisor


def schedule_from_trajectory(trajectory):
    schedule=[]
    for segment in trajectory['segments']:
        h=F(segment['duration_atomic_time']);u=[F(x) for x in segment['controls_Ha']]
        if schedule and schedule[-1][1]==u:
            schedule[-1]=(schedule[-1][0]+h,u)
        else:
            schedule.append((h,u))
    return schedule


def chebyshev_power_coefficients(values,denominator):
    """Round Chebyshev coefficients first, then change basis with integers."""
    count,rank=values.shape
    angles=(np.arange(count)+.5)*np.pi/count
    coefficients=(2/count)*np.cos(np.outer(np.arange(count),angles))@values
    coefficients[0]/=2
    rounded=[[[int(round(z.real*denominator)),int(round(z.imag*denominator))] for z in row]
             for row in coefficients]
    polynomials=[[1],[-1,2]]
    for k in range(2,count):
        p=[0]*(k+1)
        for i,a in enumerate(polynomials[-1]):
            p[i]-=2*a;p[i+1]+=4*a
        for i,a in enumerate(polynomials[-2]):
            p[i]-=a
        polynomials.append(p)
    power=[[[0,0] for _ in range(rank)] for _ in range(count)]
    for k in range(count):
        for i,weight in enumerate(polynomials[k]):
            for j,(real,imag) in enumerate(rounded[k]):
                power[i][j][0]+=weight*real;power[i][j][1]+=weight*imag
    return power


def projected_generators(data,proposal):
    states=proposal['configurations']
    v=np.asarray(proposal['vectors'],dtype=float)/proposal['denominator']
    gram=v.T@v;index={s:i for i,s in enumerate(states)}
    generators=[]
    for poly in [decode_model(data)]+[decode(c['operator'],data['modes'],4) for c in proposal['controls']]:
        d,actions=integer_action(poly,states)
        hv=np.zeros_like(v)
        for j,column in enumerate(actions):
            for state,c in column.items():
                if state in index:
                    hv[index[state]]+=(c/d)*v[j]
        generators.append(np.linalg.solve(gram,v.T@hv))
    return v,generators


def propose(data,proposal,policy,order=10,step=F(1,4),stable=False,chebyshev=False):
    start=time.monotonic();rank=len(proposal['vectors'][0])
    _,generators=projected_generators(data,proposal)
    schedule=POLICIES[policy] if isinstance(policy,str) else policy
    shift=F(round(generators[0][0,0]*10**10),10**10)
    denominator=10**14;state=np.zeros(rank,dtype=complex);state[0]=1;segments=[]
    if chebyshev:
        from scipy.linalg import expm
    for duration,control in schedule:
        remaining=duration
        a=generators[0]-float(shift)*np.eye(rank)+sum(float(u)*g for u,g in zip(control,generators[1:]))
        local_step=stable_step(a,step) if stable else step
        exponentials={}
        while remaining:
            h=min(remaining,local_step);remaining-=h
            if chebyshev:
                if h not in exponentials:
                    nodes=(1+np.cos((np.arange(order+1)+.5)*np.pi/(order+1)))/2
                    exponentials[h]=([expm(-1j*float(h)*theta*a) for theta in nodes],expm(-1j*float(h)*a))
                sample,endpoint=exponentials[h]
                values=np.asarray([e@state for e in sample])
                rounded=chebyshev_power_coefficients(values,denominator)
                state=endpoint@state
            else:
                coefficients=[state.copy()]
                for k in range(1,order+1):
                    coefficients.append((-1j*float(h)/k)*(a@coefficients[-1]))
                coefficients=np.asarray(coefficients)
                rounded=[[[int(round(z.real*denominator)),int(round(z.imag*denominator))] for z in row] for row in coefficients]
                state=coefficients.sum(axis=0)
            segments.append({'duration_atomic_time':str(h),'controls_Ha':[str(u) for u in control],'coefficients':rounded})
    result={'kind':'rational_molecular_trajectory_v1','proposal_sha256':digest(proposal),
        'coefficient_denominator':denominator,'horizon_atomic_time':str(sum(t for t,_ in schedule)),
        'phase_shift_Ha':str(shift),'segments':segments,
        'observed_spatial_orbital':min(data['modes']//2-1,data['particles']//2)}
    return result,{'policy':policy if isinstance(policy,str) else 'searched rational schedule',
        'schedule':[[str(t),[str(x) for x in u]] for t,u in schedule],
        'order':order,'maximum_step':str(step),'stable_step_rule':stable,
        'polynomial_proposal':'rounded Chebyshev interpolation' if chebyshev else 'Taylor',
        'minimum_used_step':str(min(F(s['duration_atomic_time']) for s in segments)),
        'segments':len(segments),'reduced_dimension':rank,
        'proposal_seconds':time.monotonic()-start,'floating_dynamics_is_proposal_only':True}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('fixture',type=Path);p.add_argument('proposal',type=Path)
    p.add_argument('output',type=Path);p.add_argument('--policy',choices=POLICIES,default='heldout_switch')
    p.add_argument('--schedule-from',type=Path,help='Preserve the pulse schedule of an existing trajectory')
    p.add_argument('--stable',action='store_true',help='Use order 18 and reduce step from projected generator norm')
    p.add_argument('--chebyshev',action='store_true',help='Use order-24 interpolation and exact integer polynomial basis conversion')
    a=p.parse_args()
    if a.output.exists():
        raise FileExistsError(a.output)
    policy=schedule_from_trajectory(json.loads(a.schedule_from.read_text())) if a.schedule_from else a.policy
    result,report=propose(json.loads(a.fixture.read_text()),json.loads(a.proposal.read_text()),policy,
        order=24 if a.chebyshev else (18 if a.stable else 10),stable=a.stable,chebyshev=a.chebyshev)
    a.output.write_text(json.dumps(result,separators=(',',':'))+'\n')
    a.output.with_suffix('.proposal.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report),flush=True)
