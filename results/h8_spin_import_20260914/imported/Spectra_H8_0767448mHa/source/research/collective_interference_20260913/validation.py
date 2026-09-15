"""Bounded full-sector controls, strictly separate from compact replay."""
from fractions import Fraction as F
from pathlib import Path
import json
import time

from experiments.marginal_symbolic import decode,mono,add,scale,product
from experiments.marginal_transfer_verify import apply_word
from research.collective_interference_20260913.collective import parameters,states


def full_polynomial(data,bins=1,envelope=None):
    p=parameters(data,bins);m=p['physical_modes'];B=p['count']
    if m>12:raise ValueError('Full validation capped at twelve physical modes')
    h=decode(data['active_hamiltonian'],4,4)
    na=add(*(mono(((1,i),(0,i))) for i in range(4)),mono((),-data['active_number_offset']))
    nb=add(*(mono(((1,i),(0,i))) for i in range(4,m)),mono((),-B))
    h=add(h,scale(product(na,na),F(data['charge_u'])),scale(product(nb,nb),F(data['charge_v'])),scale(product(na,nb),F(data['charge_w'])))
    gap=F(data['gap']);spread=F(data['spread']);t=p['hyb_per_orbital']
    if envelope is not None and envelope not in ('lo','hi'):raise ValueError('Invalid envelope')
    if envelope is not None:h=add(h,mono((),p['filled_energy']))
    for sign,offset,key in ((-1,4,'negative_channel'),(1,4+B,'positive_channel')):
        for j in range(B):
            mode=offset+j
            if envelope is None:e=gap+(spread*j/(B-1) if B>1 else F(0))
            else:
                g=next(g for g in p['groups'] if g['sign']==sign and g['index_first']<=j<=g['index_last'])
                e=g[envelope]
                if sign<0:h=add(h,mono((),e))
            h=add(h,mono(((1,mode),(0,mode)),sign*e))
            for i,c in enumerate(map(F,data[key])):
                if c:h=add(h,mono(((1,i),(0,mode)),t*c),mono(((1,mode),(0,i)),t*c))
    return h,p


def sparse_matrix(h,m,n):
    import numpy as np
    from scipy.sparse import coo_matrix
    basis=states(m,n);index={s:i for i,s in enumerate(basis)};entries={};work=0
    for word,c in h.items():
        for col,s in enumerate(basis):
            work+=1;target=apply_word(word,s)
            if target is not None:
                t,sign=target;key=(index[t],col);entries[key]=entries.get(key,F(0))+c*sign
    rows=[];cols=[];values=[]
    for (i,j),c in entries.items():
        if c:rows.append(i);cols.append(j);values.append(float(c))
    matrix=coo_matrix((values,(rows,cols)),shape=(len(basis),len(basis))).tocsr()
    return matrix,basis,work


def exact_sandwich_diagonals(original,bound,m,n,orientation):
    delta=add(original,scale(bound,-1)) if orientation=='lower' else add(bound,scale(original,-1))
    # With identical couplings/charging terms, only diagonal occupation terms remain.
    if any(tuple(i for f,i in w if f)!=tuple(i for f,i in w if not f) for w in delta):
        raise AssertionError('Envelope changed an off-diagonal coupling')
    minimum=None
    for s in states(m,n):
        value=F(0)
        for w,c in delta.items():
            action=apply_word(w,s)
            if action:
                target,sign=action
                if target!=s:raise AssertionError('Non-diagonal sandwich difference')
                value+=c*sign
        if value<0:raise AssertionError('Operator envelope inequality failed')
        minimum=value if minimum is None else min(minimum,value)
    return str(minimum)


def run(data,bins=1):
    import numpy as np
    from scipy.sparse.linalg import eigsh
    from research.collective_interference_20260913.collective import core_forms
    start=time.monotonic();original,p=full_polynomial(data,bins);m=p['physical_modes'];n=p['particles']
    eigenvalues={};work=0;diagonal_checks={};residuals={}
    for label in ('original','lo','hi'):
        h=original if label=='original' else full_polynomial(data,bins,label)[0]
        matrix,basis,cost=sparse_matrix(h,m,n);work+=cost
        if (matrix-matrix.T).nnz:raise AssertionError('Non-Hermitian full validation matrix')
        if matrix.shape[0]<=32:
            vals,vecs=np.linalg.eigh(matrix.toarray());ev=float(vals[0]);vector=vecs[:,0]
        else:
            vals,vecs=eigsh(matrix,k=1,which='SA',tol=1e-11,v0=np.random.default_rng(20260913).normal(size=len(basis)))
            ev=float(vals[0]);vector=vecs[:,0]
        residuals[label]=float(np.linalg.norm(matrix@vector-ev*vector))
        if residuals[label]>1e-8:raise AssertionError('Full validation eigenvector residual too large')
        eigenvalues[label]=ev
        if label!='original':
            diagonal_checks[label]=exact_sandwich_diagonals(original,h,m,n,'lower' if label=='lo' else 'upper')
            forms,_=core_forms(p,label);values=[]
            for f in forms:
                metric=np.sqrt(np.array(f['metric'],dtype=float));small=np.array(f['matrix'],dtype=float)/metric[:,None]/metric[None,:]
                values.append(float(np.linalg.eigvalsh(small)[0])+float(p['filled_energy']))
            if abs(min(values)-ev)>1e-8:raise AssertionError('Full and compressed envelope spectra disagree')
    if not eigenvalues['lo']-1e-8<=eigenvalues['original']<=eigenvalues['hi']+1e-8:raise AssertionError('Numerical spectrum escaped envelope')
    return {'physical_modes':m,'particles':n,'full_sector_dimension':len(basis),'word_state_actions':work,
            'numeric_eigenvalues':eigenvalues,'eigenvector_residuals':residuals,'exact_diagonal_sandwich_minima':diagonal_checks,
            'full_vs_compressed_envelopes_agree_1e_8':True,'wall_seconds':time.monotonic()-start,
            'scope':'Explicit finite-sector validation oracle. Numerical eigenvalues do not accept certificates.'}


if __name__=='__main__':
    import argparse
    from research.collective_interference_20260913.collective import model,propose,replay
    parser=argparse.ArgumentParser();parser.add_argument('--out',type=Path,required=True);args=parser.parse_args()
    rows=[]
    cases=[('no_dark_modes',model(1),1),('degenerate',model(2,spread='0'),1),
           ('dispersed',model(2),1),('refined',model(2),2),
           ('small_gap',model(2,gap='1/20'),1),
           ('offset_zero',{**model(2),'active_number_offset':0},1)]
    for name,data,bins in cases:
        row=run(data,bins);cert,_=propose(data,bins);accepted=replay(data,cert)
        energy=row['numeric_eigenvalues']['original']
        if not float(F(accepted['lower']))-1e-8<=energy<=float(F(accepted['upper']))+1e-8:
            raise AssertionError('Full physical ground energy escaped accepted interval')
        row.update(name=name,accepted_lower=accepted['lower'],accepted_upper=accepted['upper'])
        rows.append(row);print(name,row['full_sector_dimension'],row['numeric_eigenvalues'],flush=True)
    args.out.write_text(json.dumps({'rows':rows},indent=2)+'\n')
