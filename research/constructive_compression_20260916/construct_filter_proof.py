"""Numerical proposal only, using local terms and an explicitly charged upper."""
import argparse,json,time
from fractions import Fraction as F
from math import ceil
from pathlib import Path
import numpy as np
from research.constructive_compression_20260916.model import mpo,numerical_mpo,trace_seed,exact_seed_cert
from research.constructive_compression_20260916 import charge_mps as cm
from research.constructive_compression_20260916.export_upper import export


def from_cert(cert):
    qs=cert['bond_charges'];arrays=[];den=cert['denominator']
    for i,edges in enumerate(cert['tensors']):
        a=np.zeros((len(qs[i]),2,len(qs[i+1])))
        for l,s,r,v in edges:a[l,s,r]=v/den
        arrays.append(a)
    return arrays,[[tuple(q) for q in layer] for layer in qs]

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True)
    p.add_argument('--upper',type=Path,required=True);p.add_argument('--bond',type=int,default=16)
    p.add_argument('--steps',type=int,default=8);p.add_argument('--bits',type=int,default=32);a=p.parse_args()
    a.out.mkdir(parents=True,exist_ok=False);start=time.monotonic();rungs=4
    u=json.loads((a.upper/'upper_result.json').read_text());ell=F(ceil(F(u['upper_t'])*10**8),10**8)
    b=F(52);delta=F(9,10000);theta=np.arccosh(1+2*float(delta/(b-ell)))
    z=F(round(float(np.exp(-theta))*10**9),10**9)
    W,wq=numerical_mpo(mpo(rungs));seed=exact_seed_cert(rungs,1<<a.bits);states=[seed]
    previous=from_cert(seed);older=None;alpha=float((b+ell)/(b-ell));beta=float(-2/(b-ell))
    for i in range(1,a.steps+1):
        hp=cm.apply_mpo(*previous,W,wq)
        raw=cm.linear_combination([previous,hp] if older is None else [previous,hp,older],
          [alpha,beta] if older is None else [2*alpha,2*beta,-1])
        aa,qq,d=cm.compress(*raw,a.bond)
        data,cert=export({'arrays':[x.tolist() for x in aa],'charges':[[list(q) for q in layer] for layer in qq]},rungs,a.bits)
        states.append(cert);older,previous=previous,from_cert(cert)
        print(json.dumps({'step':i,'proposal_discarded_norm':float(np.sqrt(d['discarded_squared_norm'])),
                         'entries':sum(map(len,cert['tensors']))}),flush=True)
    proof={'kind':'hubbard_trace_chebyshev_proposal_v1','rungs':rungs,'b':str(b),'ell':str(ell),'z':str(z),
       'states':states,'upper_dependency':str(a.upper.resolve()),'proposal_bond_cap':a.bond,
       'proposal_seconds':time.monotonic()-start,'enumerated_configurations':0}
    (a.out/'proof.json').write_text(json.dumps(proof)+'\n')
