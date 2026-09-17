"""Enumerated diagnostic of local dressing on every interacting ladder bond.

No global-state fitting: all gates use the isolated t=1, U=8 bond angle.
Nonadjacent fermion gates include exact occupation-permutation signs.
"""
import argparse
from fractions import Fraction as F
import json
from pathlib import Path
import time

import numpy as np
from scipy.sparse import coo_matrix
from research.constructive_response_20260916.interacting_probe import (
    rational_rows,numerical,smallest,doublons,local_rotation,analyze,
    RATIONAL_HALF_ANGLE,
)


def permutation_sign(label,order):
    occupied=[i for i in order if label & (1<<i)]
    inversions=sum(a>b for i,a in enumerate(occupied) for b in occupied[i+1:])
    return -1 if inversions%2 else 1


def fermion_gate(labels,sites,a,b,gate):
    if not 0<=a<b<sites:
        raise ValueError('Ordered distinct sites required')
    modes=[2*a,2*a+1,2*b,2*b+1]
    order=modes+[i for i in range(2*sites) if i not in modes]
    mask=sum(1<<i for i in modes)
    lookup={v:i for i,v in enumerate(labels)}
    rr=[];cc=[];vv=[]
    for j,label in enumerate(labels):
        local=sum(((label>>mode)&1)<<k for k,mode in enumerate(modes))
        sign_in=permutation_sign(label,order)
        for target in range(16):
            x=gate[target][local]
            if not x:
                continue
            out=(label & ~mask)+sum(((target>>k)&1)<<mode for k,mode in enumerate(modes))
            sign=sign_in*permutation_sign(out,order)
            rr.append(lookup[out]);cc.append(j);vv.append(float(x)*sign)
    return coo_matrix((vv,(rr,cc)),shape=(len(labels),len(labels))).tocsr()


def circuit(rungs):
    result=[(2*i,2*i+1,'rung') for i in range(rungs)]
    for parity in (0,1):
        for i in range(parity,rungs-1,2):
            result.extend([(2*i+s,2*(i+1)+s,'horizontal_'+str(parity)) for s in (0,1)])
    return result


def run(rungs,max_nonzeros=2000000):
    start=time.monotonic()
    labels,rows=rational_rows(rungs,F(1));h=numerical(rows)
    e0,psi=smallest(h,vector=True)
    pidx=np.array([i for i,v in enumerate(labels) if doublons(v,2*rungs)==0])
    gate=local_rotation();history=[]
    for a,b,color in circuit(rungs):
        G=fermion_gate(labels,2*rungs,a,b,gate)
        h=(G.T@h@G).tocsr();h.eliminate_zeros();psi=G.T@psi
        row={'sites':[a,b],'layer':color,'matrix_nonzeros':int(h.nnz),
             'retained_weight':float(psi[pidx]@psi[pidx]),
             'ground_residual':float(np.linalg.norm(h@psi-e0*psi))}
        history.append(row)
        if h.nnz>max_nonzeros:
            return {'status':'stopped_at_declared_matrix_nonzero_cap','history':history,
                    'sites':2*rungs,'seconds':time.monotonic()-start,
                    'global_determinants_constructed':len(labels)}
    result=analyze(h,pidx,psi,e0,0.001,32)
    return {'status':'enumerated_float64_diagnostic_not_certified','sites':2*rungs,
            'U_over_t':8,'lambda':1,'global_determinants_constructed':len(labels),
            'gate_half_angle':str(RATIONAL_HALF_ANGLE),'global_state_fitting':False,
            'color_layers':3,'gate_count':len(history),'history':history,
            'response':result,'seconds':time.monotonic()-start}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--rungs',type=int,choices=[2,3,4],required=True)
    p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    if a.out.exists():
        raise FileExistsError(a.out)
    result=run(a.rungs);a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
