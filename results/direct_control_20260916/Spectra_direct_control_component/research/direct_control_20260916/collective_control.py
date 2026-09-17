"""Direct observable-only control certificate for the supplied molecular model.

Exact rational CAR/Pauli algebra and MPS scalar contractions only. A stronger,
shorter pulse is a separately declared regime, not a solution of the frozen
amplitude-.5 four-phase protocol. No trajectory or sector embedding is input.
"""
import argparse,json,time,hashlib
from collections import defaultdict
from fractions import Fraction as F
from pathlib import Path
from experiments.marginal_symbolic import decode,encode,hermitian
from experiments.marginal_hunt_car import add,mul,scale
from research.correlated_pair_20260913.mps_exact import State
from research.molecular_collective_20260913.core import digest
from .tensor_operator import product_terms


def comm(a,b):return add(mul(a,b),scale(mul(b,a),-1))


def real_pauli(poly,modes):
    """I,X,J=-iY,Z product basis. Each product has operator norm one."""
    terms=product_terms({'modes':modes,'hamiltonian':encode(poly)})
    out=defaultdict(F)
    for word,c in terms.items():
        partial={():c}
        for a,b,c0,d in word:
            choices=[('I',F(a+d,2)),('Z',F(a-d,2)),('X',F(b+c0,2)),('J',F(c0-b,2))]
            partial={w+(s,):v*k for w,v in partial.items() for s,k in choices if k}
        for w,v in partial.items():out[w]+=v
    return {w:c for w,c in out.items() if c}


def trig_interval(x,sine=False,terms=32):
    """Rational alternating Taylor enclosure after terms become decreasing."""
    x=F(x);term=x if sine else F(1);total=term
    for k in range(1,terms):
        a=2*k+(1 if sine else 0);term*=-x*x/F(a*(a-1));total+=term
    a=2*terms+(1 if sine else 0)
    nxt=term*(-x*x/F(a*(a-1)))
    if x*x >= F((a+2)*(a+1)):raise ValueError('Taylor remainder not decreasing')
    return min(total,total+nxt),max(total,total+nxt)


def atan_interval(x,terms):
    # alternating arctangent, |x|<1, x>=0
    if not 0<=x<1:raise ValueError('Arctangent domain')
    s=sum(((-1)**k*x**(2*k+1)/F(2*k+1) for k in range(terms)),F())
    nxt=(-1)**terms*x**(2*terms+1)/F(2*terms+1)
    return min(s,s+nxt),max(s,s+nxt)


def pi_interval():
    a,b=atan_interval(F(1,5),24);c,d=atan_interval(F(1,239),8)
    return 16*a-4*d,16*b-4*c


def controls(p=6,q=10):
    if any(type(x) is not int or x<0 or x%2 for x in (p,q)) or p==q:
        raise ValueError('Two distinct spatial orbitals, specified by even spin indices, are required')
    d={((1,i),(0,i)):F(c) for i,c in [(p,1),(p+1,1),(q,-1),(q+1,-1)]}
    w={((1,i),(0,j)):F(1) for i,j in [(p,q),(q,p),(p+1,q+1),(q+1,p+1)]}
    a={((1,i),(0,j)):F(c) for i,j,c in [(q,p,1),(p,q,-1),(q+1,p+1,1),(p+1,q+1,-1)]}
    return d,w,a


def construct(data,cert,target=F(-3,5),max_amplitude=F(1024),angle=F(25,8),amp_error=F(1,1000),initial_radius=F(1,2000),noise=F(1,1000),pair=(6,10)):
    start=time.monotonic()
    if max_amplitude<=0 or amp_error<0 or not 0<=initial_radius<=1 or noise<0:raise ValueError('Control/uncertainty domain')
    if type(data['modes']) is not int or data['modes']<4 or data['modes']%2:raise ValueError('Paired spin orbitals required')
    d,w,a=controls(*pair)
    if max(pair)+1>=data['modes']:raise ValueError('Control orbital outside model')
    h=decode(data['hamiltonian'],data['modes'],4)
    if not hermitian(h):raise ValueError('Hermitian drift required')
    if any(sum(2*c-1 for c,p in w if p%2==spin) for w in h for spin in (0,1)):raise ValueError('Spin populations not preserved')
    state=State(data,cert)
    if comm(w,d)!=scale(a,2) or comm(w,a)!=scale(d,2) or not hermitian(d) or not hermitian(w):raise AssertionError('Control Lie algebra')
    # Real input MPS and real skew-adjoint A imply <iA>=0 exactly.
    from experiments.marginal_hunt_car import adj
    if adj(a)!=scale(a,-1):raise AssertionError('Skew-adjoint A')
    initial=state.expectation(d)
    hd=real_pauli(comm(h,d),data['modes']);ha=real_pauli(comm(h,a),data['modes'])
    cd=sum(map(abs,hd.values()),F());cy=sum(map(abs,ha.values()),F())
    coslo,coshi=trig_interval(angle);sinlo,sinhi=trig_interval(angle,True)
    pilo,pihi=pi_interval()
    if not pihi/2<angle<pilo:raise ValueError('This closed integral formula requires pi/2 < angle < pi')
    ideal_lo=min(initial*coslo,initial*coshi);ideal_hi=max(initial*coslo,initial*coshi)
    attempts=[];amplitude=F(1,2)
    while amplitude<=max_amplitude:
        T=angle/(2*amplitude)
        # Exact integrals of |cos(2vt)| and |sin(2vt)| on this angle range.
        drift=(cd*(2-sinlo)+cy*(1-coslo))/(2*amplitude)
        robust=16*T*amp_error+4*initial_radius+4*noise
        lo=max(F(-2),ideal_lo-drift-robust);hi=min(F(2),ideal_hi+drift+robust)
        entry={'v_Ha':str(amplitude),'u_Ha':'0','T_au':str(T),'drift_allowance':str(drift),
               'robustness_allowance':str(robust),'D_interval':[str(lo),str(hi)],'D_interval_float':[float(lo),float(hi)],'target_proved':hi<=target}
        attempts.append(entry)
        if entry['target_proved']:break
        amplitude*=2
    return {'kind':'exact_collective_observable_control_v1','fixture_sha256':digest(data),'state_sha256':digest(cert),
        'initial_D':str(initial),'initial_D_float':float(initial),'initial_Y':'0','target':str(target),
        'controls':{'D':encode(d),'W':encode(w),'A_Y_divided_by_i':encode(a)},
        'HD_pauli_terms':len(hd),'HY_pauli_terms':len(ha),'HD_norm_bound_Ha':str(cd),'HY_norm_bound_Ha':str(cy),
        'angle':str(angle),'requested_max_amplitude_Ha':str(max_amplitude),'cos_interval':[str(coslo),str(coshi)],'sin_interval':[str(sinlo),str(sinhi)],
        'attempts':attempts,'accepted':attempts[-1]['target_proved'] if attempts else False,
        'uncertainties':{'per_control_pointwise_Ha':str(amp_error),'initial_trace_distance':str(initial_radius),'integrated_phase_flip_rate':str(noise)},
        'stats':state.stats,'seconds':time.monotonic()-start,
        'scope':'Direct strong-control regime on supplied molecular Hamiltonian and initial MPS. Does not solve the original H8 amplitude-.5 four-phase protocol. No laboratory actuator/preparation claim.',
        'dependencies':'Supplied rational Hamiltonian and initial MPS; their discovery is inherited and not included in response construction timing.'}


def replay(data,cert,artifact,expected_max_amplitude,expected_pair=(6,10)):
    if F(artifact['requested_max_amplitude_Ha'])!=expected_max_amplitude:raise ValueError('Declared amplitude budget mismatch')
    # The task constants and uncertainty budgets are fixed here. Saved success
    # flags, moments, commutator bounds, and endpoints never replace arithmetic.
    fresh=construct(data,cert,max_amplitude=expected_max_amplitude,pair=expected_pair)
    keys=['kind','fixture_sha256','state_sha256','initial_D','initial_Y','target','controls',
          'HD_pauli_terms','HY_pauli_terms','HD_norm_bound_Ha','HY_norm_bound_Ha','angle',
          'requested_max_amplitude_Ha','cos_interval','sin_interval','uncertainties','accepted']
    if any(fresh[k]!=artifact[k] for k in keys):raise ValueError('Collective certificate did not reproduce')
    if len(fresh['attempts'])!=len(artifact['attempts']):raise ValueError('Attempt count')
    for fresh_row,old in zip(fresh['attempts'],artifact['attempts']):
        if any(fresh_row[k]!=old[k] for k in fresh_row if k!='D_interval_float'):raise ValueError('Control or endpoint mismatch')
    return fresh


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('fixture');p.add_argument('state');p.add_argument('--max-amplitude',default='1024');p.add_argument('--output',required=True);args=p.parse_args()
    out=Path(args.output)
    if out.exists():raise FileExistsError(out)
    result=construct(json.loads(Path(args.fixture).read_text()),json.loads(Path(args.state).read_text()),max_amplitude=F(args.max_amplitude))
    out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('accepted','seconds','initial_D_float','HD_norm_bound_Ha','HY_norm_bound_Ha','HD_pauli_terms','HY_pauli_terms')}))
    print(json.dumps(result['attempts'][-1]))
