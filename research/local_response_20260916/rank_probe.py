"""Local-only numerical screening of finite-rank projector corrections."""
import argparse
import json
from pathlib import Path
from fractions import Fraction as F
from math import sqrt
import numpy as np
from research.local_response_20260916.local_exact import local_blocks


def screen(certificate, ranks=(1,2,4,6,8,12,16)):
    blocks=local_blocks(certificate['messages'])
    shifts=list(map(F,certificate['lower_shifts']))
    spectra=[]
    for bi,local in enumerate(blocks):
        values=[]
        for key,(labels,ham) in local.items():
            ev,vec=np.linalg.eigh(np.array(ham,dtype=float)-float(shifts[bi])*np.eye(len(labels)))
            for j,e in enumerate(ev):
                v=np.zeros(256);v[labels]=vec[:,j]
                values.append((float(e),v,key))
        values.sort(key=lambda x:x[0]);spectra.append(values)
    records=[]
    for rank in ranks:
        tensors=[np.stack([row[1] for row in spec[:rank]],axis=-1).reshape(16,16,rank).transpose(1,0,2)
                 for spec in spectra]
        # Contraction of three LOCAL low subspaces, no full state vector.
        k=np.einsum('abi,bcj,cdk->adjik',*tensors,optimize=True).reshape(256*rank,rank*rank)
        p=float(np.linalg.eigvalsh(k.T@k)[-1]);p=max(0,min(1,p))
        gaps=[spec[rank][0]-1e-7 for spec in spectra]
        g=min(gaps[0],gaps[2]);h=gaps[1]
        gamma=(g+h-sqrt((g-h)**2+4*g*h*p))/2
        records.append({'rank_per_patch':rank,'gap_proposals':gaps,'overlap_squared':p,
                        'correction_proposal_t':gamma,'overlap_gram_dimension':rank*rank,
                        'overlap_contraction_entries':int(k.size),
                        'global_determinants_enumerated':0,
                        'status':'numerical_proposal_only'})
    return records,spectra

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--certificate',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    if a.out.exists():raise FileExistsError(a.out)
    records,_=screen(json.loads(a.certificate.read_text()))
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(records,indent=2)+'\n')
    print(json.dumps(records,indent=2))
