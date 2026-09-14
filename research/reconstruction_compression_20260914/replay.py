"""Independent standard-library-only acceptance of expanded lower factors."""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import sys
import time
from experiments.marginal_symbolic import verify,decode

def check(fixture,certificate):
    start=time.monotonic()
    if (certificate['modes'],certificate['particles'])!=(fixture['modes'],fixture['particles']):raise ValueError('Lower sector mismatch')
    if decode(certificate['hamiltonian'],certificate['modes'],4)!=decode(fixture['hamiltonian'],fixture['modes'],4):raise ValueError('Lower Hamiltonian mismatch')
    if certificate.get('operator_degree',3)!=3:raise ValueError('Cubic certificate required')
    receipt=verify(certificate)
    receipt['replay_seconds']=time.monotonic()-start
    receipt['acceptance']='Exact CAR expansion, integer square factors, exact sector multiplier, coefficient L1 residual'
    return receipt

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('fixture');p.add_argument('certificate');p.add_argument('output');p.add_argument('--upper-receipt');a=p.parse_args()
    data=json.loads(Path(a.fixture).read_text());raw=Path(a.certificate).read_bytes();cert=json.loads(raw)
    rec=check(data,cert);rec['certificate_sha256']=hashlib.sha256(raw).hexdigest()
    if a.upper_receipt:
        upper=json.loads(Path(a.upper_receipt).read_text());U=F(upper['upper_Ha']);L=F(rec['lower'])
        if L>U:raise ValueError('Inconsistent endpoints')
        rec.update(upper_Ha=str(U),width_Ha=str(U-L),width_mHa=float(1000*(U-L)),target_1p6mHa_met=U-L<=F(1,625))
    if any(m in sys.modules for m in ('numpy','scipy','cvxpy','pyscf','quimb')):raise AssertionError('Numerical dependency in lower acceptance')
    Path(a.output).write_text(json.dumps(rec,indent=2)+'\n');print(json.dumps(rec,indent=2))
