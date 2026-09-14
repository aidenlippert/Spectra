"""Select quartic square directions from negative moment eigenvectors."""
from fractions import Fraction as F
import json
from pathlib import Path

import numpy as np

from experiments.marginal_adaptive import assemble_and_solve
from experiments.marginal_coefficient import export, gram_map, symmetric_upper
from experiments.marginal_collective import hopping_polynomial, verify_interval
from experiments.marginal_dual_exact import full_cubic_blocks, full_word_blocks


def run(iterations=6,batch=8,out_name='residual_penalty'):
    modes=10;t=F(1,5);h=hopping_polynomial(modes,t)
    groups=[[i,i+5] for i in range(5)]
    blocks=full_cubic_blocks(modes);base_count=len(blocks);candidates=full_word_blocks(modes,4)
    maps=None;history=[]
    out=Path(__file__).resolve().parents[1]/'results/marginal_cut_search'/out_name;out.mkdir(parents=True,exist_ok=True)
    for iteration in range(iterations+1):
        proposal,data=assemble_and_solve(h,modes,5,blocks,degree=4,groups=groups,residual_penalty=True)
        cert,receipt=export(h,modes,5,blocks,proposal,operator_degree=4)
        cert['variational_upper']=symmetric_upper(modes,t);receipt.update(verify_interval(cert))
        if maps is None:
            lookup={w:i for i,w in enumerate(data['rows'])}
            maps=[gram_map(b['words'],lookup) for b in candidates]
        violations=[]
        for idx,(g,b) in enumerate(zip(maps,candidates)):
            n=len(b['words']);raw=(g.T@data['dual']).reshape(n,n)
            eigenvalues,eigenvectors=np.linalg.eigh((raw+raw.T)/2)
            if eigenvalues[0]<-1e-7:
                violations.append((float(eigenvalues[0]),idx,eigenvectors[:,0]))
        violations.sort(key=lambda x:x[0])
        receipt.update({'iteration':iteration,'added_rank_one_blocks':len(blocks)-base_count,
                        'negative_candidate_blocks':len(violations),'worst_eigenvalue':violations[0][0] if violations else 0})
        (out/f'iteration{iteration}.json').write_text(json.dumps(cert)+'\n')
        (out/f'iteration{iteration}_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
        history.append(receipt);(out/'history.json').write_text(json.dumps(history,indent=2)+'\n')
        print(json.dumps({k:receipt[k] for k in ('iteration','lower_float','width_float','worst_eigenvalue','solve_seconds')}),flush=True)
        if iteration==iterations or not violations:break
        for value,idx,vector in violations[:batch]:
            candidate=candidates[idx]
            blocks.append({'name':f'cut{iteration}:{idx}','words':candidate['words'],
                           'basis_transform':vector[:,None]})
    return history


if __name__=='__main__':run()
