"""Exact limitation of the single-pulse global commutator-sum envelope.

One explicit number/spin-compatible determinant is a norm witness. No sector
basis is generated. This is not a statement of physical unreachability.
"""
from fractions import Fraction as F
from collections import defaultdict
import argparse,json
from pathlib import Path
from experiments.marginal_symbolic import decode,encode,hermitian
from .collective_control import controls,comm,trig_interval
from research.intervention_reduction_20260916.exact import sqrt_up
from research.correlated_pair_20260913.mps_exact import State
from research.molecular_collective_20260913.core import digest


def all_angle_floor(initial,cd_lower,cy_lower,max_amplitude):
    """Lower bound on the upper endpoint of this envelope, for every duration.

    Positive initial D and zero initial Y. Ignore nonnegative uncertainty
    allowances. Split theta=2|v|T into [0,pi/2], [pi/2,pi], and [pi,infinity).
    This is a limitation of the norm-sum proof rule, not actual reachability.
    """
    if initial<0 or min(cd_lower,cy_lower)<0 or max_amplitude<=0:raise ValueError('Envelope domain')
    a=cd_lower/(2*max_amplitude);b=cy_lower/(2*max_amplitude)
    middle=2*a+b-sqrt_up((initial-b)**2+a*a,scale=1<<80)
    late=2*a+2*b-initial
    return min(F(),middle,late),[F(),middle,late]


def selected_action(poly,label):
    out=defaultdict(F)
    for word,c in poly.items():
        state=label;sign=1
        for create,p in reversed(word):
            if ((state>>p)&1)==create:sign=0;break
            if (state & ((1<<p)-1)).bit_count()%2:sign=-sign
            state^=1<<p
        if sign:out[state]+=c*sign
    return {k:v for k,v in out.items() if v}


def calculate(data,cert,control_receipt,amplitude=F(1,2)):
    if data['modes']!=16 or cert['spin_counts']!=[4,4]:raise ValueError('H8 obstruction domain')
    h=decode(data['hamiltonian'],data['modes'],4);d,w,a=controls()
    if not hermitian(h):raise ValueError('Hermitian drift required')
    if control_receipt['fixture_sha256']!=digest(data) or control_receipt['state_sha256']!=digest(cert):raise ValueError('Obstruction input binding')
    if control_receipt['controls']!={'D':encode(d),'W':encode(w),'A_Y_divided_by_i':encode(a)}:raise ValueError('Obstruction control binding')
    state=State(data,cert);initial=state.expectation(d)
    # State accepts real integer tensors only. A is real and skew-adjoint,
    # so <iA>=0 exactly, rather than an unverified receipt assumption.
    if control_receipt['initial_Y']!='0' or F(control_receipt['initial_D'])!=initial:raise ValueError('Obstruction initial moments')
    label=(1<<data['particles'])-1
    if [(label & sum(1<<p for p in range(s,data['modes'],2))).bit_count() for s in (0,1)]!=[4,4]:raise ValueError('Witness sector')
    images=[selected_action(comm(h,x),label) for x in (d,a)]
    norms2=[sum((c*c for c in image.values()),F()) for image in images]
    lower=[max(F(),sqrt_up(x,scale=1<<80)-F(1,1<<80)) for x in norms2]
    theta=F(control_receipt['angle']);cl,ch=trig_interval(theta);sl,sh=trig_interval(theta,True)
    ideal_lo=min(initial*cl,initial*ch)
    drift_min=(lower[0]*(2-sh)+lower[1]*(1-ch))/(2*amplitude)
    T=theta/(2*amplitude);u=control_receipt['uncertainties']
    robustness=16*T*F(u['per_control_pointwise_Ha'])+4*F(u['initial_trace_distance'])+4*F(u['integrated_phase_flip_rate'])
    endpoint_floor=ideal_lo+drift_min+robustness
    every_angle,ranges=all_angle_floor(initial,lower[0],lower[1],amplitude)
    return {'kind':'exact_commutator_sum_envelope_obstruction_v1','witness_determinant':label,
       'fixture_sha256':digest(data),'state_sha256':digest(cert),'initial_moments_recomputed':True,
       'input_determinants_used':1,'full_sector_enumeration':False,'action_image_supports':[len(x) for x in images],
       'commutator_norm_squared_lower':[str(x) for x in norms2],'commutator_norm_lower':[str(x) for x in lower],
       'amplitude_upper_Ha':str(amplitude),'angle':str(theta),'unavoidable_drift_allowance_at_max_amplitude':str(drift_min),
       'best_envelope_upper_endpoint_is_at_least':str(endpoint_floor),'endpoint_floor_float':float(endpoint_floor),
       'target':control_receipt['target'],'family_cannot_certify':endpoint_floor>F(control_receipt['target']),
       'all_durations_upper_endpoint_floor':str(every_angle),'angle_range_floors':[str(x) for x in ranges],
       'all_durations_cannot_certify':every_angle>F(control_receipt['target']),
       'scope':'Single constant W pulse with u=0, global commutator-sum upper-envelope rule, |v| <= declared maximum. Both fixed-angle and all-duration floors are reported. Even exact operator norms cannot make this rule pass. Not physical unreachability, not the original four-phase family, and not a many-body complexity theorem.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('fixture');p.add_argument('state');p.add_argument('control');p.add_argument('output');a=p.parse_args()
    rec=calculate(json.loads(Path(a.fixture).read_text()),json.loads(Path(a.state).read_text()),json.loads(Path(a.control).read_text()))
    out=Path(a.output)
    if out.exists():raise FileExistsError(out)
    out.write_text(json.dumps(rec,indent=2)+'\n');print(json.dumps(rec))
