"""Numerical operator construction from exact finite Hubbard coefficients."""
import numpy as np
from .exact import validate_model,model_hash as identity,source,integer_hamiltonian

def compile_hubbard(spec):
    h=integer_hamiltonian(spec);cores=[]
    for i,entries in enumerate(h['cores']):
        ops={}
        for u,v,p,q,c in entries:
            if (u,v) not in ops:ops[u,v]=np.zeros((4,4))
            # Exactly one GLOBAL denominator: placed on the first core.
            ops[u,v][p,q]+=c/h['den'] if i==0 else c
        cores.append(ops)
    return dict(cores=cores,charges=h['charges'],max_operator_bond=max(map(len,h['charges'])))
