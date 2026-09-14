"""Exact balanced-spin proof, lifted to the full even-N molecular sector.

The numerical proposer is not imported. The existing CAR checker verifies
an auxiliary two-body Hamiltonian on the whole N sector. It agrees with the
exact SU(2) average of the input H on M_S=0. Every even-N spin multiplet has
an M_S=0 member; the coefficient norm pays the original H's spin defect.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import sys
import time
from experiments.marginal_symbolic import (decode,encode,mono,add,scale,product,
                                           hermitian,verified_residual)
from research.certificate_scaling.spin_twirl import twirl,ladder

def alpha_shift(m,n):
    return add(mono((),-F(n,2)),*(mono(((1,i),(0,i))) for i in range(0,m,2)))

def setup(data):
    m,n=data['modes'],data['particles']
    if type(m) is not int or type(n) is not int or m<2 or m%2 or n%2 or not 0<=n<=m:
        raise ValueError('Even modes and even particle number required')
    h=decode(data['hamiltonian'],m,4)
    if not hermitian(h) or any(sum(2*c-1 for c,i in w) for w in h):
        raise ValueError('Hermitian number-conserving Hamiltonian required')
    hs=twirl(h)
    if ladder(hs,True) or ladder(hs,False):raise AssertionError('SU(2) projection failed')
    delta=sum(abs(c) for c in add(h,scale(hs,-1)).values())
    return h,hs,delta

def check(data,cert):
    start=time.monotonic();h,hs,delta=setup(data);m,n=data['modes'],data['particles']
    if cert.get('kind')!='balanced_spin_sos_v1':raise ValueError('Unknown certificate kind')
    if (cert.get('modes'),cert.get('particles'))!=(m,n):raise ValueError('Sector mismatch')
    if decode(cert['hamiltonian'],m,4)!=h:raise ValueError('Input Hamiltonian mismatch')
    y=decode(cert['alpha_multiplier'],m,2)
    if not hermitian(y) or any(sum(2*c-1 for c,i in w) or sum(2*c-1 for c,i in w if i%2==0) for w in y):
        raise ValueError('Spin multiplier must preserve both spin counts and be Hermitian')
    core=cert['core']
    if core.get('operator_degree')!=3 or (core['modes'],core['particles'])!=(m,n):raise ValueError('Core sector or degree mismatch')
    expected=add(hs,scale(product(alpha_shift(m,n),y),-1))
    if decode(core['hamiltonian'],m,4)!=expected:raise ValueError('Auxiliary Hamiltonian identity mismatch')
    residual,receipt=verified_residual(core)
    lower=F(receipt['lower'])-delta
    receipt.update(kind=cert['kind'],lower=str(lower),lower_float=float(lower),
        original_H_spin_defect_Ha=str(delta),
        balanced_spin_counts=[n//2,n//2],valid_on='Entire fixed-N sector',
        lifting='Exact SU(2) projection, every even-N multiplet meets M_S=0, original-H coefficient-norm perturbation',
        many_body_states_enumerated=0,replay_seconds=time.monotonic()-start)
    return receipt

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('fixture');p.add_argument('certificate');p.add_argument('output');p.add_argument('--upper');a=p.parse_args()
    data=json.loads(Path(a.fixture).read_text());raw=Path(a.certificate).read_bytes();cert=json.loads(raw)
    rec=check(data,cert);rec['certificate_sha256']=hashlib.sha256(raw).hexdigest()
    if a.upper:
        U=F(json.loads(Path(a.upper).read_text())['upper_Ha']);L=F(rec['lower'])
        if U<L:raise ValueError('Inconsistent endpoints')
        rec.update(upper_Ha=str(U),width_Ha=str(U-L),width_mHa=float(1000*(U-L)),target_1p6mHa_met=U-L<=F(1,625))
    if any(k in sys.modules for k in ('numpy','scipy','cvxpy','pyscf','quimb')):raise AssertionError('Numerical acceptance dependency')
    with Path(a.output).open('x') as f:json.dump(rec,f,indent=2)
    print(json.dumps(rec,indent=2))
