"""The full accepting-L1 objective on the unchanged Gram family.

This numerical conic formulation includes every Gram entry and every residual
coefficient. Its output still requires the original rational certificate checker.
"""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys
import time
import numpy as np
from scipy import sparse
from research.acceptance_channels_20260915.campaign import ROOT
from research.acceptance_channels_20260915.accepted_fit import residual_lift
from research.acceptance_channels_20260915.scs_control import standard_problem,pack,unpack
from research.gpu_acceleration_20260915.solve import Operator
from research.gpu_acceleration_20260915.kernel_bench import cpu_project
from research.sector_quotient_20260914.search import export


def accepted_problem(op,ids,K,weights):
    """Dual: min h.y, F'.y=e0, A'.y>=0, y=K'.d, |d|<=w."""
    base,cone=standard_problem(op,ids)
    nr,nf=K.shape[1],K.shape[0]
    if len(weights)!=nf or np.any(np.asarray(weights)<=0):
        raise ValueError('Positive residual weights and matching lift required')
    zero=sparse.csc_matrix((len(ids),nf))
    constraints=[sparse.hstack((base['A'][:len(ids)],zero)),
        sparse.hstack((sparse.eye(nr),-K.T)),
        sparse.hstack((sparse.csc_matrix((nf,nr)),sparse.eye(nf))),
        sparse.hstack((sparse.csc_matrix((nf,nr)),-sparse.eye(nf))),
        sparse.hstack((base['A'][len(ids):],sparse.csc_matrix((base['A'].shape[0]-len(ids),nf))))]
    data={'A':sparse.vstack(constraints,format='csc'),
        'b':np.r_[base['b'][:len(ids)],np.zeros(nr),weights,weights,np.zeros(base['A'].shape[0]-len(ids))],
        'c':np.r_[op.rhs,np.zeros(nf)]}
    return data,{'z':len(ids)+nr,'l':2*nf,'s':cone['s']}


def run(case,source,tag,seconds):
    import scs
    started=time.monotonic();out=case/tag;out.mkdir(exist_ok=False)
    sys.path.insert(0,str(ROOT/'.venv-interacting-libs'))
    from research.nvidia_followup_20260915.sparse_quotient import SparseQuotient
    op=Operator(case);fixture=case/'fixture.json'
    if hashlib.sha256(fixture.read_bytes()).hexdigest()!=op.meta['fixture_sha256']:
        raise ValueError('Hamiltonian changed')
    op.meta['fixture']=str(fixture)
    quo=SparseQuotient(op);ids=np.r_[0,quo.pivots]
    raw=np.load(source);oldQ=[raw[f'Q_{i}'] for i in range(len(op.Q))];oldx=quo.recover(oldQ)
    K,lift_details=residual_lift(op.prepared,op.scale)
    weights=np.load(op.prepared/'weights.npy');nr,nf=len(op.rhs),len(weights)
    data,cone=accepted_problem(op,ids,K,weights)
    moment=-raw['y']
    T=sparse.load_npz(op.prepared/'twirl.npz');selected=np.load(op.prepared/'selected.npy')
    d=T[selected].T@(op.scale*moment)
    warm_x=np.r_[moment,np.clip(d,-weights,weights)]
    r=op.rhs-op.A(oldQ)-op.free@oldx;full=K@r
    warm_y=np.r_[-oldx[ids],-r,np.maximum(-full,0),np.maximum(full,0),*(pack(Q) for Q in oldQ)]
    warm_s=data['b']-data['A']@warm_x
    setup_seconds=time.monotonic()-started
    print(json.dumps({'stage':'full_accepted_L1_SCS','A_shape':list(data['A'].shape),
        'A_nonzeros':data['A'].nnz,'setup_seconds':setup_seconds,'version':scs.__version__}),flush=True)
    solver=scs.SCS(data,cone,time_limit_secs=seconds,max_iters=20000,
        eps_abs=1e-8,eps_rel=1e-8,acceleration_lookback=10,verbose=True)
    result=solver.solve(warm_start=True,x=warm_x,y=warm_y,s=warm_s)
    if any(not np.isfinite(result[k]).all() for k in ('x','y','s')):
        raise ValueError('No finite accepted-objective iterate')
    offset=len(ids)+nr+2*nf;Q=[]
    for old in op.Q:
        n=len(old);count=n*(n+1)//2
        Q.append(unpack(result['y'][offset:offset+count],n));offset+=count
    Q=cpu_project(Q,'evd')
    x=np.zeros(op.free.shape[1]);x[ids]=-result['y'][:len(ids)]
    def score(Q,x):
        eta=float(weights@abs(K@(op.rhs-op.A(Q)-op.free@x)))
        return float(x[0]-eta),eta
    proposals=[('source',oldQ,oldx),('accepted_objective',Q,x),('least_squares_refit',Q,quo.recover(Q))]
    scored=[(score(q,z),name,q,z) for name,q,z in proposals]
    (lower,eta),chosen,Q,x=max(scored,key=lambda p:p[0][0])
    exported=export(op,Q,x,out/'export')
    np.savez_compressed(out/'checkpoint.npz',x=x,y=-result['x'][:nr],**{f'Q_{i}':q for i,q in enumerate(Q)})
    U=float(Fraction(op.meta['upper_Ha']));delta=float(Fraction(op.meta['original_H_spin_defect_Ha']))
    record={'kind':'full_accepted_L1_SCS','same_Gram_family':True,'all_Gram_variables_free':True,
        'source':str(source),'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
        'setup_seconds':setup_seconds,'total_seconds':time.monotonic()-started,
        'SCS_info':{k:v.item() if isinstance(v,np.generic) else v for k,v in result['info'].items()},
        'selected_proposal':chosen,'proposal_predicted_lower_Ha':{name:sc[0] for sc,name,_,_ in scored},
        'predicted_width_mHa':1000*(U-lower+delta),'predicted_residual_Ha':eta,
        'exact_family_obstruction':False,'exact_acceptance_required':True,'lift':lift_details,'export':exported}
    (out/'discovery.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('case',type=Path);p.add_argument('source',type=Path)
    p.add_argument('tag');p.add_argument('--seconds',type=float,default=180)
    a=p.parse_args();run(a.case.resolve(),a.source.resolve(),a.tag,a.seconds)
