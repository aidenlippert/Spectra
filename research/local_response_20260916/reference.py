"""Explicitly charged global enumeration for a physical upper/reference ONLY."""
import argparse
from fractions import Fraction as F
import json
from pathlib import Path
from itertools import combinations


def check_upper(payload):
    if payload.get('kind')!='enumerated_open_2x4_hubbard_upper_v1':
        raise ValueError('Wrong fixed model')
    choices=list(combinations(range(8),4))
    labels=sorted(sum(1<<(2*i) for i in up)+sum(1<<(2*i+1) for i in down)
                  for up in choices for down in choices)
    if payload.get('labels')!=labels:raise ValueError('Wrong fixed charge ordering')
    z=payload['integer_vector']
    if len(z)!=len(labels) or any(type(x) is not int for x in z):
        raise ValueError('Integer vector does not match the sector')
    norm=sum(x*x for x in z)
    if not norm:raise ValueError('Zero vector')
    index={x:i for i,x in enumerate(labels)}
    edges=[(2*i,2*i+1) for i in range(4)]
    edges += [(2*i+s,2*(i+1)+s) for i in range(3) for s in (0,1)]
    energy=0;entries=0
    for col,state in enumerate(labels):
        doublons=sum(((state>>(2*i))&3)==3 for i in range(8))
        energy+=8*doublons*z[col]*z[col]
        entries+=int(bool(doublons))
        for a,b in edges:
            for spin in (0,1):
                for dst,src in ((2*a+spin,2*b+spin),(2*b+spin,2*a+spin)):
                    if state&(1<<src) and not state&(1<<dst):
                        lo,hi=sorted((dst,src))
                        between=((1<<hi)-1)^((1<<(lo+1))-1)
                        sign=-1 if (state&between).bit_count()%2 else 1
                        target=state^(1<<src)^(1<<dst)
                        energy-=sign*z[index[target]]*z[col]
                        entries+=1
    return {'status':'accepted_exact_enumerated_variational_upper',
            'energy_upper_over_t':str(F(energy,norm)),
            'global_determinants_enumerated':len(labels),
            'matrix_entries_evaluated':entries,
            'scope':'Physical upper only; never used to discover local lower or projector correction'}


if __name__=='__main__':
    import numpy as np
    from research.constructive_response_20260916.interacting_probe import rational_rows,numerical,smallest
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    if a.out.exists():raise FileExistsError(a.out)
    a.out.mkdir(parents=True)
    labels,rows=rational_rows(4,F(1));h=numerical(rows);e,psi=smallest(h,vector=True)
    payload={'kind':'enumerated_open_2x4_hubbard_upper_v1','labels':labels,
             'integer_vector':[int(round(x*10**10)) for x in psi]}
    receipt=check_upper(payload)
    receipt['numerical_reference_energy']=e
    receipt['numerical_eigen_residual']=float(np.linalg.norm(h@psi-e*psi))
    (a.out/'upper.json').write_text(json.dumps(payload,separators=(',',':'))+'\n')
    (a.out/'construction.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,indent=2))
