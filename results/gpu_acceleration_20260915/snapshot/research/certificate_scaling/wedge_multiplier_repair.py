"""Fixed-factor sparse LP for number-ideal repair with a wedge residual bound."""
import argparse,json,sys,time
from fractions import Fraction as F
from itertools import combinations
from math import comb
from pathlib import Path
import numpy as np
from scipy import sparse
from scipy.optimize import linprog
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from experiments.marginal_symbolic import add,decode,encode,expand_squares,multiplier_basis,number_shift,product,scale
from research.certificate_scaling.wedge_residual_bound import replay
from research.certificate_scaling.adaptive_block_discovery import partition
from research.certificate_scaling.direct_sparse_discovery import sparse_columns

def repair(cert,seconds=60):
    start=time.monotonic();before=replay(cert);m,n=cert['modes'],cert['particles']
    if cert.get('operator_degree',3)!=3:raise ValueError('Repair supports cubic factors')
    h=decode(cert['hamiltonian'],m,4);sq,_=expand_squares(cert['blocks'],cert['denominator'],m)
    target=add(h,scale(sq,-1))
    _,signature,_=partition(h,m,'quadratic',True);zero=signature(())
    basis=[p for p in multiplier_basis(m,max_body=2) if all(signature(w)==zero for w in p)]
    cols=[product(number_shift(m,n),q) for q in basis]
    words=set(target)|{w for p in cols for w in p}
    rows=sorted((w for w in words if len(w)//2<=n and tuple(i for c,i in w if c)<=tuple(i for c,i in w if not c)),key=lambda w:(len(w),w))
    lookup={w:i for i,w in enumerate(rows)}
    signs=np.array([(-1)**((len(w)//2)*(len(w)//2-1)//2) for w in rows])
    a=sparse.diags(signs)@sparse_columns([{w:c for w,c in p.items() if w in lookup} for p in cols],lookup)
    a=a.tocsr();rhs=np.array([float(target.get(w,0)) for w in rows])*signs
    levels=list(range(1,min(3,n)+1));diagonals=[(k,indices) for k in levels for indices in combinations(range(m),k)]
    didx={v:i for i,v in enumerate(diagonals)}
    off=[i for i,w in enumerate(rows) if tuple(p for c,p in w if c)!=tuple(p for c,p in w if not c)]
    select_r=[];select_c=[];inc_r=[];inc_c=[]
    for i,(k,indices) in enumerate(diagonals):
        w=tuple((1,p) for p in indices)+tuple((0,p) for p in indices)
        if w in lookup:select_r.append(i);select_c.append(lookup[w])
    for j,ri in enumerate(off):
        w=rows[ri];k=len(w)//2
        for indices in (tuple(p for c,p in w if c),tuple(p for c,p in w if not c)):
            inc_r.append(didx[k,indices]);inc_c.append(j)
    select=sparse.csr_matrix((np.ones(len(select_r)),(select_r,select_c)),shape=(len(diagonals),len(rows)))
    incidence=sparse.csr_matrix((np.ones(len(inc_r)),(inc_r,inc_c)),shape=(len(diagonals),len(off)))
    level_map=sparse.csr_matrix((np.ones(len(diagonals)),(range(len(diagonals)),[levels.index(k) for k,_ in diagonals])),shape=(len(diagonals),len(levels)))
    ao=a[off];ident=sparse.eye(len(off));zeros=sparse.csr_matrix((len(off),len(levels)))
    inequalities=sparse.vstack([sparse.hstack([-ao,-ident,zeros]),sparse.hstack([ao,-ident,zeros]),sparse.hstack([select@a,incidence,level_map])],format='csc')
    bounds_rhs=np.r_[-rhs[off],rhs[off],select@rhs]
    constant=a[lookup[()]].toarray().ravel() if () in lookup else np.zeros(len(basis))
    objective=np.r_[constant,np.zeros(len(off)),[-comb(n,k) for k in levels]]
    built=time.monotonic()
    solved=linprog(objective,A_ub=inequalities,b_ub=bounds_rhs,bounds=[(None,None)]*len(basis)+[(0,None)]*len(off)+[(None,None)]*len(levels),method='highs',options={'time_limit':seconds})
    solved_at=time.monotonic();out=cert;after=before
    if solved.x is not None and np.all(np.isfinite(solved.x)):
        candidate=dict(cert);candidate['number_multiplier']=encode(add(*(scale(p,F(round(float(v)*10**12),10**12)) for p,v in zip(basis,solved.x[:len(basis)]))))
        candidate_bound=replay(candidate)
        if F(candidate_bound['lower'])>=F(before['lower']):out=candidate;after=candidate_bound
    return out,{'before':before,'after':after,'retained_improvement':out is not cert,'LP_status':int(solved.status),'LP_message':solved.message,'LP_iterations':solved.nit,'ideal_variables':len(basis),'offdiagonal_variables':len(off),'wedge_rows':len(diagonals),'LP_nonzeros':inequalities.nnz,'build_seconds':built-start,'solve_seconds':solved_at-built,'wall_seconds':time.monotonic()-start,'factors_changed':False,'objective':'sum of exact fixed-N wedge Gershgorin lower bounds; constant residual treated exactly'}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--input',type=Path,required=True);p.add_argument('--outputdir',type=Path,required=True);p.add_argument('--seconds',type=float,default=60);args=p.parse_args()
    cert,rec=repair(json.loads(args.input.read_text()),args.seconds);args.outputdir.mkdir(parents=True,exist_ok=True)
    (args.outputdir/'certificate.json').write_text(json.dumps(cert,separators=(',',':'))+'\n')
    (args.outputdir/'receipt.json').write_text(json.dumps(rec,indent=2)+'\n');print(json.dumps(rec),flush=True)
