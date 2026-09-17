"""Independent numerical algorithm for the unchanged equality-constrained SDP.

The PSD packing follows https://www.cvxgrp.org/scs/api/cones.html.
This is not a new certificate checker or a full coefficient-L1 epigraph solver.
"""
import argparse
from fractions import Fraction
import hashlib,json,sys,time
from pathlib import Path
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import splu
from research.acceptance_channels_20260915.campaign import ROOT
from research.gpu_acceleration_20260915.solve import Operator
from research.gpu_acceleration_20260915.kernel_bench import cpu_project
from research.sector_quotient_20260914.search import export


def coordinates(n):
    i=np.array([i for j in range(n) for i in range(j,n)],dtype=int)
    j=np.array([j for j in range(n) for i in range(j,n)],dtype=int)
    return i,j,np.where(i==j,1.,np.sqrt(2.))


def pack(Q):
    i,j,w=coordinates(len(Q));return Q[i,j]*w


def unpack(v,n):
    i,j,w=coordinates(n);Q=np.zeros((n,n));Q[i,j]=v/w;Q[j,i]=v/w;return Q


def standard_problem(op,ids):
    rows=[op.free[:,ids].T.tocsc()]
    for M,Q in zip(op.M,op.Q):
        i,j,w=coordinates(len(Q))
        packed=M[:,i*len(Q)+j]@sparse.diags(w)
        rows.append(-packed.T.tocsc())
    A=sparse.vstack(rows,format='csc')
    target=np.zeros(A.shape[0]);target[0]=1
    return {'A':A,'b':target,'c':op.rhs.copy()},{'z':len(ids),'s':[len(Q) for Q in op.Q]}


def run(case,source,tag,seconds):
    import scs
    started=time.monotonic();out=case/tag;out.mkdir(exist_ok=False)
    sys.path.insert(0,str(ROOT/'.venv-interacting-libs'))
    from research.nvidia_followup_20260915.sparse_quotient import SparseQuotient
    op=Operator(case);fixture=case/'fixture.json'
    if hashlib.sha256(fixture.read_bytes()).hexdigest()!=op.meta['fixture_sha256']:raise ValueError('Hamiltonian changed')
    op.meta['fixture']=str(fixture)
    quo=SparseQuotient(op);ids=np.r_[0,quo.pivots]
    raw=np.load(source);oldQ=[raw[f'Q_{i}'] for i in range(len(op.Q))]
    oldx=quo.recover(oldQ)
    data,cone=standard_problem(op,ids)
    x0=-raw['y'];y0=np.r_[-oldx[ids],*(pack(Q) for Q in oldQ)]
    s0=data['b']-data['A']@x0
    build_seconds=time.monotonic()-started
    print(json.dumps({'stage':'SCS','A_shape':list(data['A'].shape),'A_nonzeros':data['A'].nnz,
        'setup_seconds':build_seconds,'version':scs.__version__}),flush=True)
    solver=scs.SCS(data,cone,time_limit_secs=seconds,max_iters=20000,
        eps_abs=1e-8,eps_rel=1e-8,acceleration_lookback=10,verbose=True)
    result=solver.solve(warm_start=True,x=x0,y=y0,s=s0)
    if any(not np.isfinite(result[k]).all() for k in ('x','y','s')):
        raise ValueError('No finite SCS iterate')
    offset=len(ids);Q=[]
    for old in op.Q:
        n=len(old);count=n*(n+1)//2
        Q.append(unpack(result['y'][offset:offset+count],n));offset+=count
    Q=cpu_project(Q,'evd')
    x=quo.recover(Q)
    T=sparse.load_npz(op.prepared/'twirl.npz');selected=np.load(op.prepared/'selected.npy')
    lift=T[:,selected];inverse=splu(T[selected][:,selected].tocsc())
    weights=np.load(op.prepared/'weights.npy')
    def score(Q,x):
        r=(op.rhs-op.A(Q)-op.free@x)/op.scale
        eta=float(weights@abs(lift@inverse.solve(r)))
        return float(x[0]-eta),eta
    old_score,old_eta=score(oldQ,oldx);new_score,eta=score(Q,x)
    kept_new=new_score>old_score
    if not kept_new:Q,x,eta=oldQ,oldx,old_eta
    exported=export(op,Q,x,out/'export')
    np.savez_compressed(out/'checkpoint.npz',x=x,y=-result['x'],**{f'Q_{i}':q for i,q in enumerate(Q)})
    U=float(Fraction(op.meta['upper_Ha']));delta=float(Fraction(op.meta['original_H_spin_defect_Ha']))
    receipt={'kind':'unchanged_family_SCS_control','same_Gram_family':True,
        'source':str(source),'setup_seconds':build_seconds,'total_seconds':time.monotonic()-started,
        'SCS_info':{k:v.item() if isinstance(v,np.generic) else v for k,v in result['info'].items()},
        'new_iterate_selected':kept_new,'predicted_width_mHa':1000*(U-x[0]+eta+delta),
        'predicted_residual_Ha':eta,'exact_family_obstruction':False,
        'exact_acceptance_required':True,'export':exported}
    (out/'discovery.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('case',type=Path);p.add_argument('source',type=Path)
    p.add_argument('tag');p.add_argument('--seconds',type=float,default=180)
    a=p.parse_args();run(a.case.resolve(),a.source.resolve(),a.tag,a.seconds)
