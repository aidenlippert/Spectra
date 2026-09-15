"""Independent numerical CAR matrix assembly for small benchmark references only."""
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import argparse
import json
import time

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh

from research.side_routes_20260913.chain_campaign import case


def reference(task):
    name,m = task
    start = time.monotonic()
    model = case(name,m); n = model['particles']
    if m > 16: raise ValueError('Enumeration reference is limited to 16 modes')
    states = [sum(1<<i for i in indices) for indices in combinations(range(m),n)]
    position = {s:i for i,s in enumerate(states)}
    rows,cols,values = [],[],[]
    t,v,h = [list(map(F,model[key])) for key in ('hopping','interaction','fields')]
    for col,mask in enumerate(states):
        bits = [(mask>>i)&1 for i in range(m)]
        diagonal = sum(h[i]*bits[i] for i in range(m))+sum(v[i]*bits[i]*bits[i+1] for i in range(m-1))
        rows.append(col); cols.append(col); values.append(float(diagonal))
        for bond in range(m-1):
            for create,annihilate in ((bond,bond+1),(bond+1,bond)):
                state = mask; sign = 1
                for is_create,site in ((False,annihilate),(True,create)):
                    occupied = bool(state&(1<<site))
                    if occupied == is_create:
                        sign = 0; break
                    if (state&((1<<site)-1)).bit_count()%2: sign = -sign
                    state ^= 1<<site
                if sign:
                    rows.append(position[state]); cols.append(col); values.append(float(-t[bond]*sign))
    matrix = coo_matrix((values,(rows,cols)),shape=(len(states),len(states))).tocsr()
    if (matrix-matrix.T).nnz: raise AssertionError('Independent CAR Hamiltonian is not symmetric')
    energy,vector = eigsh(matrix,k=1,which='SA',tol=1e-12,v0=np.ones(len(states))/np.sqrt(len(states)))
    energy = float(energy[0]); vector = vector[:,0]
    residual = float(np.linalg.norm(matrix@vector-energy*vector))
    if residual > 1e-9: raise AssertionError('Numerical reference did not satisfy residual gate')
    return {'case':name,'modes':m,'model':model,'dimension':len(states),'matrix_nonzeros':matrix.nnz,
            'energy_float':energy,'residual_norm':residual,'wall_seconds':time.monotonic()-start,
            'scope':'Numerical Lanczos reference, independently assembled using CAR. Enumerates the small fixed-N basis; not a formal lower certificate or part of discovery.'}


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('--out',type=Path,required=True); a=p.parse_args()
    started=time.monotonic()
    with ProcessPoolExecutor(max_workers=3) as pool:
        results=list(pool.map(reference,[(name,m) for name in ('repulsive','bond_disorder','site_disorder') for m in (8,12,16)]))
    receipt={'references':results,'wall_seconds':time.monotonic()-started,'workers':3}
    a.out.write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'references':len(results),'wall_seconds':receipt['wall_seconds']}))
