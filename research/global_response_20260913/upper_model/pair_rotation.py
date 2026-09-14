"""Compact variational upper candidate from pair rotations of an HF determinant.

The ansatz is a product of occupied/virtual one-particle rotations.  Its
finite support is generated from the six rotation parameters, never from a
reference many-body eigenvector.  The exact checker uses CAR word actions.
"""
from pathlib import Path
from fractions import Fraction
import json, math, time, sys
import numpy as np
from scipy.optimize import minimize
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from experiments.marginal_symbolic import decode
from research.certificate_scaling.streaming_reference_upper import compile_term, upper

ROOT=Path(__file__).resolve().parents[3]
OUT=ROOT/'results/global_response_20260913/upper_model'

def terms(data):
    return [compile_term(w,float(c)) for w,c in decode(data['hamiltonian'],data['modes'],4).items() if compile_term(w,float(c))]

def ansatz(modes, particles, angles):
    # Interleaved spin orbitals: rotate occupied spatial orbitals 0..n-1 into
    # virtual partners n..2n-1, independently for both spins.
    n=modes//2; ne=particles//2
    out={(1<<(2*ne))-1:1.0};
    for spin in (0,1):
        for i,a in enumerate(angles[spin*ne:(spin+1)*ne]):
            occ=2*i+spin; virt=2*(ne+i)+spin
            nxt={}
            for s,x in out.items():
                if not (s>>occ)&1: raise ValueError('HF occupancy convention failed')
                nxt[s]=nxt.get(s,0)+x*math.cos(a)
                t=s^(1<<occ)^(1<<virt)
                # creation replacement sign from ordered CAR operators
                sign=(-1)**((s & ((1<<virt)-1)).bit_count() + (s & ((1<<occ)-1)).bit_count())
                nxt[t]=nxt.get(t,0)+x*math.sin(a)*sign
            out=nxt
    q=math.sqrt(sum(x*x for x in out.values()))
    return {s:x/q for s,x in out.items() if abs(x/q)>1e-14}

def action(v,ts,modes):
    r={}
    for s,x in v.items():
        for req,occ,flip,par,c in ts:
            if s&req==occ:
                t=s^flip; r[t]=r.get(t,0)+x*c*(-1 if (s&par).bit_count()%2 else 1)
    return r

def energy(v,ts,modes):
    hv=action(v,ts,modes)
    return sum(x*hv.get(s,0) for s,x in v.items())

def wick_expectation(word, R, modes):
    """Gaussian Wick contraction for a CAR word, with rho[p,q]=<c†p c_q>."""
    if not word: return 1.0
    if len(word)%2: return 0.0
    def pair(a,b):
        ca,ia=a; cb,ib=b
        if ca and cb: return 0.0
        if ca and not cb: return R[ib,ia]
        if not ca and cb: return (1.0 if ia==ib else 0.0)-R[ia,ib]
        return 0.0
    def rec(seq):
        if not seq: return 1.0
        a=seq[0]; total=0.0
        for j in range(1,len(seq)):
            total += (-1)**(j-1)*pair(a,seq[j])*rec(seq[1:j]+seq[j+1:])
        return total
    return rec(tuple(word))

def wick_energy(data, angles):
    modes=data['modes']; n=modes//2; ne=data['particles']//2
    # columns of occupied rotated orbitals; pair rotation uses tan(angle).
    C=np.zeros((modes,data['particles']))
    col=0
    for spin in (0,1):
        for i,a in enumerate(angles[:ne] if spin==0 else angles[ne:]):
            o=2*i+spin; q=2*(ne+i)+spin; t=math.tan(a); z=math.sqrt(1+t*t)
            C[o,col]=1/z; C[q,col]=t/z; col+=1
    R=C@C.T
    total=0.0; count=0
    for w,c in decode(data['hamiltonian'],modes,4).items():
        total += float(c)*wick_expectation(w,R,modes); count+=1
    return total,count

def run_case(name):
    p=ROOT/'results/certificate_scaling/active_space_ladder'/name/'fixture.json'
    data=json.loads(p.read_text()); modes=data['modes']; particles=data['particles']; ne=particles//2
    ts=terms(data); rng=np.random.default_rng(20260913+len(name)); best=None; start=time.monotonic()
    for seed in range(3):
        x=rng.normal(0,.08,2*ne)
        r=minimize(lambda a: wick_energy(data,a)[0],x,method='BFGS',options={'maxiter':40,'gtol':1e-7})
        if best is None or r.fun<best.fun: best=r
    v=ansatz(modes,particles,best.x); den=10**12
    ints={s:int(round(a*den)) for s,a in v.items() if int(round(a*den))}
    witness={'states':sorted(ints),'amplitudes':[ints[s] for s in sorted(ints)]}
    exact,_=upper(data,witness)
    hf=wick_energy(data,np.zeros(2*ne))[0]
    we,terms_count=wick_energy(data,best.x)
    return {'case':name,'sector_dimension':data['sector_dimension'],'pair_count':2*ne,'support':len(v),'support_fraction':len(v)/data['sector_dimension'],'parameters':2*ne,'float_energy':best.fun,'wick_energy':we,'wick_terms':terms_count,'wick_vs_support_energy_error':we-energy(v,ts,modes),'hf_energy':hf,'gain_mHa':1000*(hf-we),'exact_upper_Ha':str(exact),'optimizer_success':bool(best.success),'seconds':time.monotonic()-start,'witness':witness,'source':'HF pair rotations; no prior amplitudes'}

if __name__=='__main__':
    rows=[run_case(x) for x in ('h6','h8')]
    OUT.mkdir(parents=True,exist_ok=True); (OUT/'pair_rotation_results.json').write_text(json.dumps(rows,indent=2)+'\n')
    print(json.dumps([{k:r[k] for k in r if k!='witness'} for r in rows],indent=2))
