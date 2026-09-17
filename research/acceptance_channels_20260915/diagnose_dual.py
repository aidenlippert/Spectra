"""Numerical separating-direction diagnostic; never an accepted obstruction."""
import argparse
from fractions import Fraction
import json
from pathlib import Path
import time
import numpy as np
from scipy import linalg, sparse
from scipy.sparse.linalg import LinearOperator, cg, splu
from research.gpu_acceleration_20260915.solve import Operator
from research.gpu_acceleration_20260915.kernel_bench import cpu_project


def run(case, source, gap, output):
    start=time.monotonic()
    op=Operator(case)
    data=np.load(source)
    W=[data[f'W_{i}'] for i in range(len(op.Q))]
    Q=cpu_project(W,'evd'); z=data['z']
    U=float(Fraction(op.meta['upper_Ha'])); b=U-gap/1000
    f0=op.free[:,0]; F=op.free[:,1:]
    G=(op.G-f0@f0.T).tocsc()
    pre=LinearOperator(G.shape,matvec=splu(G+sparse.eye(G.shape[0],format='csc')*1e-11).solve)
    rhs=op.rhs-b*f0.toarray().ravel()-op.A([2*q-w for q,w in zip(Q,W)])-F@z
    direction,status=cg(G,rhs,M=pre,rtol=1e-12,atol=1e-15,maxiter=150)
    normalization=float((f0.T@direction)[0])
    if abs(normalization)<1e-14 or not np.isfinite(direction).all():
        raise ValueError('No stable nonzero separating-direction normalization')
    y=direction/normalization
    eigenvalues=[float(linalg.eigvalsh(A,subset_by_index=(0,0),check_finite=False)[0]) for A in op.AT(y)]
    T=sparse.load_npz(op.prepared/'twirl.npz');selected=np.load(op.prepared/'selected.npy')
    moments=T[selected].T@(op.scale*y)/np.load(op.prepared/'weights.npy')
    energy=float(op.rhs@y)
    record={'kind':'unverified_feasibility_displacement_dual','source':str(source),
        'seconds':time.monotonic()-start,'raw_normalization':normalization,'cg_status':int(status),
        'dual_energy_Ha':energy,'unverified_interval_floor_mHa':1000*(U-energy),
        'maximum_ideal_defect':float(np.max(abs(F.T@y))),
        'minimum_Gram_eigenvalue':min(eigenvalues),'maximum_moment_absolute_value':float(np.max(abs(moments))),
        'exactly_repaired':False,'family_obstruction_proved':False}
    with output.open('x') as stream:json.dump(record,stream,indent=2);stream.write('\n')
    np.savez_compressed(output.with_suffix('.npz'),y=-y)
    print(json.dumps(record),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('case',type=Path);p.add_argument('source',type=Path)
    p.add_argument('output',type=Path);p.add_argument('--gap-mha',type=float,default=1.)
    a=p.parse_args();run(a.case.resolve(),a.source.resolve(),a.gap_mha,a.output.resolve())
