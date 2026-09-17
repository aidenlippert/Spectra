"""Compact orbit proofs plus local squares for unequal hopping perturbations.

Numerical diagonalization proposes only the upper witness. Replay is standard
library only and evaluates that integer vector against the certificate's H.
"""
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import time

from experiments.marginal_symbolic import decode, encode, add
from experiments.marginal_orbit_certificate import verify
from experiments.marginal_transfer_verify import apply_word

ROOT=Path(__file__).resolve().parents[1]


def hopping_perturbation(modes,epsilon):
    if type(modes) is not int or modes<4 or modes%2:raise ValueError('Even modes >=4 required')
    half=modes//2;delta={};squares=[];shift=F(0)
    for i in range(half):
        e=F(epsilon)*(i-F(half-1,2))
        if not e:continue
        left,right=i,i+half;sign=1 if e>0 else -1;shift+=abs(e)
        delta[((1,left),(0,right))]=-e;delta[((1,right),(0,left))]=-e
        for p in ({((0,left),):F(1),((0,right),):F(-sign)},
                  {((1,left),):F(1),((1,right),):F(sign)}):
            squares.append({'polynomial':encode(p),'weight':str(abs(e)/2)})
    return delta,squares,shift


def integer_upper(h,modes,particles):
    import numpy as np
    states=[s for s in range(1<<modes) if s.bit_count()==particles]
    index={s:i for i,s in enumerate(states)};matrix=np.zeros((len(states),len(states)))
    for word,c in h.items():
        for j,state in enumerate(states):
            result=apply_word(word,state)
            if result:
                destination,sign=result;matrix[index[destination],j]+=float(c)*sign
    if not np.allclose(matrix,matrix.T,atol=1e-12,rtol=0):raise ValueError('Non-Hermitian Hamiltonian')
    values,vectors=np.linalg.eigh(matrix)
    amplitudes=[int(round(float(x)*10**10)) for x in vectors[:,0]]
    return {'amplitudes':amplitudes,'numerical_energy':float(values[0])}


def replay(certificate):
    lower_receipt=verify(certificate)
    modes,particles=certificate['modes'],certificate['particles']
    states=[s for s in range(1<<modes) if s.bit_count()==particles]
    index={s:i for i,s in enumerate(states)}
    witness=certificate.get('independent_upper')
    if type(witness) is not dict:raise ValueError('Missing upper witness')
    amps=witness.get('amplitudes')
    if type(amps) is not list or len(amps)!=len(states) or any(type(a) is not int for a in amps) or not any(amps):
        raise ValueError('Expected nonzero integer amplitudes in ascending fixed-N bitstring order')
    h=decode(certificate['hamiltonian'],modes,4);norm=sum(a*a for a in amps);energy=F(0)
    for word,c in h.items():
        for state,a in zip(states,amps):
            if not a:continue
            result=apply_word(word,state)
            if result:
                destination,sign=result;energy+=c*a*amps[index[destination]]*sign
    upper=energy/norm;width=upper-F(lower_receipt['lower'])
    if width<0:raise ValueError('Inconsistent interval')
    return dict(lower_receipt,upper=str(upper),upper_float=float(upper),norm=str(norm),
                width=str(width),width_float=float(width),hamiltonian_bound=True)


def run(epsilon=F(1,1000),source=None):
    source=ROOT/'results/marginal_distill_direct/compact_certificate.json' if source is None else Path(source)
    started=time.monotonic();raw=source.read_bytes();base=json.loads(raw)
    modes,particles=base['modes'],base['particles']
    h=decode(base['hamiltonian'],modes,4);delta,squares,shift=hopping_perturbation(modes,epsilon)
    new_h=add(h,delta);upper=integer_upper(new_h,modes,particles)
    baseline=dict(base,hamiltonian=encode(new_h),independent_upper=upper)
    baseline.pop('variational_upper',None)
    hybrid=dict(baseline,b=str(F(base['b'])-shift),
                direct_squares=list(base.get('direct_squares',[]))+squares)
    baseline_receipt=replay(baseline);hybrid_receipt=replay(hybrid)
    receipt={'epsilon':str(epsilon),'source':str(source),'source_sha256':hashlib.sha256(raw).hexdigest(),
             'perturbation':'Centered unequal pair hoppings epsilon*(i-(half-1)/2)',
             'baseline':baseline_receipt,'hybrid':hybrid_receipt,'added_direct_squares':len(squares),
             'analytic_lower_shift':str(-shift),'seconds':time.monotonic()-started,
             'scope':'Analytic perturbation transfer; square directions are not reoptimized for the new H.'}
    out=ROOT/'results/marginal_symmetry_transfer'/('hopping_'+str(epsilon).replace('/','_'))
    out.mkdir(parents=True,exist_ok=True)
    for name,item in [('baseline_certificate',baseline),('certificate',hybrid),('receipt',receipt)]:
        (out/f'{name}.json').write_text(json.dumps(item,indent=2)+'\n')
    print(json.dumps(receipt),flush=True)
    return receipt


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--epsilon',default='1/1000')
    parser.add_argument('--source',type=Path);parser.add_argument('--verify',type=Path)
    args=parser.parse_args()
    if args.verify:print(json.dumps(replay(json.loads(args.verify.read_text())),indent=2))
    else:run(F(args.epsilon),args.source)
