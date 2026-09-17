"""Charge the circuit action without expanding the transformed Hamiltonian.

Still enumerates a global determinant vector space; this only removes the
transformed-matrix storage failure of layered_probe, not the many-body cost.
"""
import argparse
from fractions import Fraction as F
import json
from pathlib import Path
import time
from math import sqrt
import numpy as np
from scipy.sparse.linalg import LinearOperator,eigsh
from research.constructive_response_20260916.interacting_probe import (
    rational_rows,numerical,smallest,local_rotation,doublons,opnorm_small_right,
)
from research.constructive_response_20260916.layered_probe import circuit,fermion_gate


def setup(rungs):
    labels,rows=rational_rows(rungs,F(1));h=numerical(rows)
    gates=[fermion_gate(labels,2*rungs,a,b,local_rotation()) for a,b,_ in circuit(rungs)]
    def action(x):
        y=x
        for g in reversed(gates): y=g@y
        y=h@y
        for g in gates: y=g.T@y
        return y
    return labels,h,gates,action


def run(rungs):
    start=time.monotonic();labels,h,gates,action=setup(rungs);n=len(labels)
    e0,psi=smallest(h,vector=True)
    for g in gates: psi=g.T@psi
    P=np.array([i for i,x in enumerate(labels) if doublons(x,2*rungs)==0]);r=len(P)
    keep=np.zeros(n,dtype=bool);keep[P]=True;Q=np.flatnonzero(~keep);nq=len(Q)
    def eliminated(x):
        if x.ndim==1:
            y=np.zeros(n);y[Q]=x
        else:
            y=np.zeros((n,x.shape[1]));y[Q]=x
        return action(y)[Q]
    C=LinearOperator((nq,nq),matvec=eliminated,matmat=eliminated,dtype=float)
    v0=np.sin(np.arange(nq)+0.731)
    cmin=float(eigsh(C,k=1,which='SA',v0=v0,tol=2e-11,return_eigenvectors=False)[0])
    cmax=float(eigsh(C,k=1,which='LA',v0=v0,tol=2e-11,return_eigenvectors=False)[0])
    embed=np.zeros((n,r));embed[P,np.arange(r)]=1
    hp=action(embed);B=hp[Q];A=hp[P]
    e=e0-.001;delta=cmin-e
    diag=np.array([8*doublons(labels[i],2*rungs)-e for i in Q])
    def ashift(x): return eliminated(x)-e*x
    X=np.zeros_like(B);history=[];first=None
    for step in range(33):
        R=B-ashift(X)
        if step in (0,1,2,4,8,16,32):
            norm=opnorm_small_right(R);bound=norm*norm/delta
            history.append({'iteration':step,'residual_norm':norm,'gap_bound':bound})
            if first is None and bound<=.001: first=step
        if step<32: X+=R/diag[:,None]
    norm=float(psi[P]@psi[P]);res=float(np.linalg.norm(action(psi)-e0*psi))
    return {'status':'enumerated_float64_diagnostic_not_certified','sites':2*rungs,
        'global_determinants_enumerated':n,'retained_dimension':r,
        'base_H_nonzeros':int(h.nnz),'gate_nonzeros_total':sum(int(g.nnz) for g in gates),
        'stored_transformed_H_entries':0,'response_entries':int(X.size),
        'matrix_free_does_not_mean_determinant_free':True,
        'preconditioner':'8 times bare doublon count minus target energy',
        'gate_count':len(gates),'ground_energy':e0,'ground_residual':res,
        'retained_weight':norm,'coupling_norm':opnorm_small_right(B),
        'eliminated_gap':cmin-e0,'condition_at_target':(cmax-e)/delta,
        'first_sampled_iteration_below_0p001':first,
        'residual_history':history,'seconds':time.monotonic()-start}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--rungs',type=int,choices=[2,3,4],required=True)
    p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    if a.out.exists(): raise FileExistsError(a.out)
    result=run(a.rungs);a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
