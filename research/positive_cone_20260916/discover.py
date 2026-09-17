"""Cold matrix-space candidate solve, then exact positive-cone acceptance.

Krylov iteration retains all 4,900 amplitude coordinates. No saved full ground
state or earlier interval endpoint is used to discover the candidate.
"""
import argparse
from fractions import Fraction as F
import json
from math import floor
from pathlib import Path
import time
import numpy as np
from scipy.linalg import eigh
from scipy.sparse.linalg import LinearOperator,eigsh
from research.positive_cone_20260916.cone import one_spin_model,verify


def construct(out):
    if out.exists():raise FileExistsError(out)
    out.mkdir(parents=True)
    start=time.monotonic()
    labels,k,edges,eta=one_spin_model();n=len(labels)
    K=np.array(k,dtype=float)
    potential=np.array([[8*(4-(a&b).bit_count()) for b in labels] for a in labels],dtype=float)
    calls=0
    def action(v):
        nonlocal calls
        calls+=1;C=np.reshape(v,(n,n))
        return (K@C+C@K+potential*C).reshape(-1)
    operator=LinearOperator((n*n,n*n),matvec=action,rmatvec=action,dtype=float)
    values,vectors=eigsh(operator,k=1,which='SA',v0=np.eye(n).reshape(-1),tol=1e-13,maxiter=10000,ncv=32)
    C=vectors[:,0].reshape(n,n);C=(C+C.T)/2
    if np.trace(C)<0:C=-C
    numerical_residual=float(np.linalg.norm(action(C.reshape(-1))-values[0]*C.reshape(-1)))
    min_eig_C=float(np.linalg.eigvalsh(C)[0])
    if min_eig_C<=0:raise RuntimeError('Candidate ground matrix not strictly positive')
    attempts=[]
    # Rounding and all scalar endpoints remain untrusted until integer replay.
    for scale in (10**10,10**12,10**14):
        z=np.rint(C*scale).astype(np.int64)
        z=(z+z.T)//2
        zr=z.astype(float)/scale
        Lzr=K@zr+zr@K+potential*zr
        proposed=float(eigh(Lzr,zr,eigvals_only=True,subset_by_index=[0,0])[0])
        payload={'kind':'half_filled_bipartite_hubbard_positive_cone_v1',
            'rungs':4,'U':'8','t':'1','target_spin_populations':[4,4],
            'one_spin_labels':labels,'integer_C':[[int(x) for x in row] for row in z]}
        for margin in (1e-8,1e-7,1e-6,1e-5):
            payload['lower_over_t']=str(F(floor((proposed-margin)*10**10),10**10))
            try:
                accepted=verify(payload)
                attempts.append({'amplitude_scale':scale,'lower_margin':margin,'accepted':True,
                                 'width_float':accepted['width_float']})
                if accepted['meets_target']:break
            except ValueError as error:
                attempts.append({'amplitude_scale':scale,'lower_margin':margin,'accepted':False,'reason':str(error)})
        else:continue
        if accepted['meets_target']:break
    else:raise RuntimeError('No exact interval at the frozen target')
    elapsed=time.monotonic()-start
    receipt={'acceptance':accepted,'discovery_method':'Cold full-amplitude matrix-space Krylov from identity',
             'global_state_teacher_used':False,'map_applications':calls,
             'numerical_energy':float(values[0]),'numerical_residual_frobenius':numerical_residual,
             'numerical_min_eigenvalue_of_C':min_eig_C,
             'elapsed_internal_seconds':elapsed,'repair_attempts':attempts,
             'not_a_scalability_breakthrough':True}
    (out/'certificate.json').write_text(json.dumps(payload,separators=(',',':'))+'\n')
    (out/'construction.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,indent=2),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True)
    construct(p.parse_args().out)
