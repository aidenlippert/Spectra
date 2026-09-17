"""Cold four-site block response proposal; no full eigenvector/eigenvalue solve.

Uses enumerated 30-by-6 response columns and six-dimensional Ritz solves. All
floating choices are checked by the separate rational accepting checker.
"""
import argparse
from fractions import Fraction as F
import json
from math import ceil,floor,sqrt
from pathlib import Path
import time
import numpy as np
from scipy.linalg import eigh
from research.constructive_response_20260916.metric_certificate import square_model,check_square


def propose():
    start=time.monotonic();labels,H,P,Q=square_model()
    h=np.array(H,dtype=float);C=h[np.ix_(Q,Q)];B=h[np.ix_(Q,P)];A=h[np.ix_(P,P)]
    e=-1.0;X=np.zeros_like(B);history=[]
    for outer in range(1,9):
        D=C-e*np.eye(len(Q));diag=np.diag(D)
        for _ in range(8):
            X += (B-D@X)/diag[:,None]
        X=np.round(X*1e8)/1e8
        M=np.eye(len(P))+X.T@X
        G=A-B.T@X-X.T@B+X.T@C@X
        vals,vecs=eigh((G+G.T)/2,M)
        e=floor((float(vals[0])-1e-9)*1e7)/1e7
        v=vecs[:,0]
        R=B-(C-e*np.eye(len(Q)))@X
        relative2=float(eigh(R.T@R,M,eigvals_only=True)[-1])
        absolute2=float(np.linalg.eigvalsh(R.T@R)[-1])
        rho=F(ceil(sqrt(max(0,relative2))*1e7)+1,10**7)
        absrho=F(ceil(sqrt(max(0,absolute2))*1e7)+1,10**7)
        payload={'kind':'square_hubbard_metric_response_v1',
                 'reference_energy_over_t':str(F(round(e*1e7),10**7)),
                 'rho_over_t':str(rho),'absolute_rho_over_t':str(absrho),
                 'response':[[str(F(round(x*1e8),10**8)) for x in row] for row in X],
                 'retained_trial':[str(F(round(x*1e8),10**8)) for x in v]}
        accepted=check_square(payload)
        width=F(accepted['width_over_t'])
        history.append({'outer_iteration':outer,'jacobi_steps_cumulative':8*outer,
                        'reference_energy_over_t':e,'accepted_width_over_t':str(width)})
        if width<=F(1,1000):
            return payload,{'seconds':time.monotonic()-start,'history':history,
                'discovery':'30-by-6 enumerated block Jacobi, six-dimensional Ritz solves',
                'global_eigensolver_calls':0,'global_determinants_enumerated':len(labels),
                'successful_interval':accepted}
    raise RuntimeError('Bounded response construction did not meet target')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    if a.out.exists():
        raise FileExistsError(a.out)
    a.out.mkdir(parents=True)
    cert,record=propose()
    (a.out/'certificate.json').write_text(json.dumps(cert,indent=2)+'\n')
    (a.out/'construction.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(record,indent=2))
