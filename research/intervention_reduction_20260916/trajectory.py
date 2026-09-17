"""Direct exact residual acceptance for a molecular state trajectory.

Checks rational polynomial segments in an integer subspace using CAR actions.
The residual integral is computed jointly: cancellations between basis columns,
Hamiltonian terms and derivatives are retained. No floating integration enters
acceptance, and every segment jump is paid.
"""
import argparse
from fractions import Fraction as F
from math import comb,lcm
import json
from pathlib import Path
import resource
import sys
import time

from research.intervention_reduction_20260916.exact import (
    digest, exact_moments, rational, sqrt_up,
)


def validate_coefficients(coefficients, rank):
    if not coefficients or len(coefficients)>33:
        raise ValueError('Polynomial orders zero through 32 supported')
    for vector in coefficients:
        if len(vector)!=rank or any(not isinstance(z,list) or len(z)!=2 or any(type(x) is not int for x in z) for z in vector):
            raise ValueError('Complex integer coefficient dimensions')


def linear_combination(rows, coefficients):
    return {s:(sum(x*z[0] for x,z in zip(row,coefficients)),
               sum(x*z[1] for x,z in zip(row,coefficients))) for s,row in rows.items()}


def vector_norm_square(amplitudes):
    return sum(a*a+b*b for a,b in amplitudes.values())


def polynomial_norm_range(gram,coefficients,coefficient_den,initial_norm):
    """Exact Bernstein enclosure of ||Vp(theta)||^2/G00, theta in [0,1]."""
    degree=2*(len(coefficients)-1)
    applied=[[[sum(g*z[part] for g,z in zip(row,c)) for part in (0,1)]
              for row in gram] for c in coefficients]
    power=[0]*(degree+1)
    for i,c in enumerate(coefficients):
        for j in range(i,len(coefficients)):
            dot=sum(a*x+b*y for (a,b),(x,y) in zip(c,applied[j]))
            power[i+j]+=(1 if i==j else 2)*dot
    scale=coefficient_den**2*initial_norm
    bernstein=[sum((F(power[i]*comb(k,i),comb(degree,i)*scale) for i in range(k+1)),F(0))
               for k in range(degree+1)]
    return min(bernstein),max(bernstein)


def mixed_action(moments, rows, controls, shift):
    weights=[F(1)]+controls
    dens=moments['action_denominators']
    common=lcm(*(d*w.denominator for d,w in zip(dens,weights)),shift.denominator)
    scales=[int(common*w/d) for d,w in zip(dens,weights)]
    actions=moments['integer_actions']
    rank=moments['rank']
    labels=set(rows).union(*(set(a) for a in actions))
    zero=[0]*rank
    result={}
    for label in labels:
        row=[sum(scale*a.get(label,zero)[j] for scale,a in zip(scales,actions))
             -int(common*shift)*rows.get(label,zero)[j] for j in range(rank)]
        if any(row):
            result[label]=row
    return common,result


def residual_integral(rows, action, action_den, coefficients, duration, coefficient_den, initial_norm):
    """Integral over theta in [0,1] of ||iVp' - duration*HVp||^2 / initial norm."""
    p=len(coefficients)-1
    physical=[linear_combination(rows,c) for c in coefficients]
    applied=[linear_combination(action,c) for c in coefficients]
    labels=set(rows)|set(action)
    derivative_scale=action_den*duration.denominator
    residual=[]
    for k in range(p+1):
        vector=[]
        for label in labels:
            ar,ai=applied[k].get(label,(0,0))
            dr,di=physical[k+1].get(label,(0,0)) if k<p else (0,0)
            vector.append((-derivative_scale*(k+1)*di-duration.numerator*ar,
                           derivative_scale*(k+1)*dr-duration.numerator*ai))
        residual.append(vector)
    integral_den=lcm(*range(1,2*p+2))
    numerator=0
    for i in range(p+1):
        for j in range(i,p+1):
            dot=sum(a*c+b*d for (a,b),(c,d) in zip(residual[i],residual[j]))
            numerator+=(1 if i==j else 2)*dot*(integral_den//(i+j+1))
    if numerator<0:
        raise AssertionError('Negative exact integrated squared residual')
    return F(numerator, integral_den*derivative_scale**2*coefficient_den**2*initial_norm)


def check_trajectory(data, proposal, trajectory, *, _kernel=None):
    started=time.monotonic()
    proposal_id=digest(proposal) if _kernel is None else _kernel.proposal_sha256
    if trajectory.get('kind')!='rational_molecular_trajectory_v1' or trajectory.get('proposal_sha256')!=proposal_id:
        raise ValueError('Trajectory binding')
    den=trajectory['coefficient_denominator']
    if type(den) is not int or den<=0:
        raise ValueError('Positive integer trajectory denominator required')
    horizon=rational(trajectory['horizon_atomic_time'])
    if horizon<=0:
        raise ValueError('Positive horizon required')
    shift=rational(trajectory['phase_shift_Ha'])
    if _kernel is None:
        moments=exact_moments(data,proposal)
        rows=dict(zip(proposal['configurations'],proposal['vectors']))
        initial_norm=sum(row[0]**2 for row in rows.values())
        modes=data['modes'];vector_denominator=proposal['denominator'];fixture_id=digest(data)
    else:
        if ((data is not None and _kernel.fixture_sha256!=digest(data))
                or (proposal is not None and _kernel.proposal_sha256!=digest(proposal))
                or _kernel.orbital!=trajectory['observed_spatial_orbital']):
            raise ValueError('Compiled kernel binding')
        moments=_kernel.moments;initial_norm=_kernel.metric[0][0]
        modes=_kernel.modes;vector_denominator=_kernel.vector_denominator;fixture_id=_kernel.fixture_sha256
    rank=moments['rank']
    integer_gram=[[int(g*vector_denominator**2) for g in row] for row in moments['metric']]
    previous=[[den if j==0 else 0,0] for j in range(rank)]
    total_time=F(0); error=F(0);jump_total=F(0);residual_total=F(0);cache={};segments=[]
    minimum_norm_square=None;maximum_norm_square=None
    if not trajectory['segments']:
        raise ValueError('At least one segment required')
    for segment in trajectory['segments']:
        duration=rational(segment['duration_atomic_time'])
        if duration<=0:
            raise ValueError('Positive segment duration required')
        controls=[rational(x) for x in segment['controls_Ha']]
        if len(controls)!=len(moments['amplitudes']) or any(abs(u)>a for u,a in zip(controls,moments['amplitudes'])):
            raise ValueError('Control outside declared amplitude box')
        coefficients=segment['coefficients'];validate_coefficients(coefficients,rank)
        difference=[[a-c,b-d] for (a,b),(c,d) in zip(coefficients[0],previous)]
        jump_integer=vector_norm_square(linear_combination(rows,difference)) if _kernel is None else _kernel.norm_square(difference)
        jump_sq=F(jump_integer,den**2*initial_norm)
        jump=sqrt_up(jump_sq)
        key=tuple(controls)
        if _kernel is None:
            if key not in cache:
                cache[key]=mixed_action(moments,rows,controls,shift)
            action_den,action=cache[key]
            integral=residual_integral(rows,action,action_den,coefficients,duration,den,initial_norm)
        else:
            integral=_kernel.residual_integral(controls,shift,coefficients,duration,den,initial_norm)
        norm_range=polynomial_norm_range(integer_gram,coefficients,den,initial_norm)
        minimum_norm_square=norm_range[0] if minimum_norm_square is None else min(minimum_norm_square,norm_range[0])
        maximum_norm_square=norm_range[1] if maximum_norm_square is None else max(maximum_norm_square,norm_range[1])
        local=sqrt_up(integral)
        error+=jump+local;jump_total+=jump;residual_total+=local;total_time+=duration
        previous=[[sum(c[j][component] for c in coefficients) for component in (0,1)] for j in range(rank)]
        segments.append({'integrated_residual_squared':str(integral),'residual_bound':str(local),'jump_bound':str(jump),
            'norm_square_enclosure':[str(x) for x in norm_range]})
    if total_time!=horizon:
        raise ValueError('Segment durations do not match horizon')
    endpoint=linear_combination(rows,previous) if _kernel is None else None
    endpoint_norm=vector_norm_square(endpoint) if _kernel is None else _kernel.norm_square(previous)
    if not endpoint_norm:
        raise ValueError('Zero endpoint state cannot define normalized observable')
    norm_ratio=F(endpoint_norm,den**2*initial_norm)
    norm_up=sqrt_up(norm_ratio);norm_lo=norm_up if norm_up**2==norm_ratio else norm_up-F(1,10**12)
    norm_correction=max(abs(1-norm_lo),abs(1-norm_up))
    normalized_error=error+norm_correction
    uniform_error=None
    if minimum_norm_square>0:
        min_norm=sqrt_up(minimum_norm_square)
        if min_norm**2!=minimum_norm_square:
            min_norm-=F(1,10**12)
        max_norm=sqrt_up(maximum_norm_square)
        uniform_error=error+max(abs(1-min_norm),abs(1-max_norm))
    orbital=trajectory['observed_spatial_orbital']
    if type(orbital) is not int or not 0<=2*orbital+1<modes:
        raise ValueError('Observed spatial orbital outside model')
    occupation=lambda s: ((s>>(2*orbital))&1)+((s>>(2*orbital+1))&1)
    if _kernel is None:
        expected=F(sum(occupation(s)*(a*a+b*b) for s,(a,b) in endpoint.items()),endpoint_norm)
        initial=F(sum(occupation(s)*row[0]**2 for s,row in rows.items()),initial_norm)
    else:
        expected=F(_kernel.observable_numerator(previous),endpoint_norm)
        initial=F(_kernel.observable[0][0],initial_norm)
    margin=2*normalized_error
    lower=max(F(0),expected-margin);upper=min(F(2),expected+margin)
    forbidden=[s for s in ('numpy','scipy','quimb','pyscf','numba') if s in sys.modules]
    if forbidden:
        raise AssertionError(('Numerical accepting dependency',forbidden))
    return {'status':'accepted_molecular_trajectory','fixture_sha256':fixture_id,'proposal_sha256':proposal_id,
        'trajectory_sha256':digest(trajectory),'horizon_atomic_time':str(horizon),'segments':segments,
        'phase_aligned_state_error_bound':str(error),'normalized_state_error_bound':str(normalized_error),
        'normalized_state_error_float':float(normalized_error),'target_met':normalized_error<=F(1,200),
        'uniform_in_time_state_error_bound':None if uniform_error is None else str(uniform_error),
        'uniform_in_time_state_error_float':None if uniform_error is None else float(uniform_error),
        'uniform_in_time_target_met':uniform_error is not None and uniform_error<=F(1,200),
        'uniform_in_time_scope':'normalized polynomial prediction on each segment, both sides of every jump, phase aligned',
        'jump_total':str(jump_total),'residual_total':str(residual_total),'normalization_correction':str(norm_correction),
        'observed_spatial_orbital':orbital,'initial_population':str(initial),'predicted_population':str(expected),
        'population_interval':[str(lower),str(upper)],'population_interval_float':[float(lower),float(upper)],
        'population_change_interval':[str(lower-initial),str(upper-initial)],
        'population_change_interval_float':[float(lower-initial),float(upper-initial)],
        'initial_state':'normalized first exact proposal column; not asserted to be the exact ground state',
        'control_scope':'only the specified piecewise-constant protocol','physical_model_error_included':False,
        'reduced_dimension':rank,**moments['counts'],'exact_replay_seconds':time.monotonic()-started,
        'peak_RSS_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'accepting_numerical_imports':forbidden}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('fixture',type=Path);p.add_argument('proposal',type=Path)
    p.add_argument('trajectory',type=Path);p.add_argument('output',type=Path);a=p.parse_args()
    if a.output.exists():
        raise FileExistsError(a.output)
    result=check_trajectory(*(json.loads(path.read_text()) for path in (a.fixture,a.proposal,a.trajectory)))
    a.output.write_text(json.dumps(result,separators=(',',':'))+'\n')
    print(json.dumps({k:result[k] for k in ('target_met','reduced_dimension','normalized_state_error_float',
        'population_change_interval_float','exact_replay_seconds','peak_RSS_bytes')}),flush=True)
