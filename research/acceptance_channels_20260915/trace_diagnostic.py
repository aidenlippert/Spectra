"""Numerical singlet-trace feasibility diagnostic, never an exact dual proof."""
import argparse
from fractions import Fraction
import json
from pathlib import Path
import time
import sys
import numpy as np
from scipy import linalg,sparse
from research.acceptance_channels_20260915.accepted_fit import residual_lift
from research.gpu_acceleration_20260915.solve import Operator
from research.interacting_scaling_20260915.singlet_trace import magnetic_trace,singlet_dimension
from research.acceptance_channels_20260915.campaign import ROOT


def run(case,checkpoint,output):
    started=time.monotonic();op=Operator(case)
    m,n=op.meta['modes'],op.meta['particles'];dimension=singlet_dimension(m,n)
    raw=np.zeros(len(op.meta['rows']))
    for i,word in enumerate(op.meta['rows']):
        w=tuple(map(tuple,word));creators=tuple(k for c,k in w if c)
        if creators!=tuple(k for c,k in w if not c):continue
        p={w:Fraction(1)}
        raw[i]=float((magnetic_trace(p,m,n//2,n//2)-magnetic_trace(p,m,n//2+1,n//2-1))/dimension)
    K,details=residual_lift(op.prepared,op.scale)
    weights=np.load(op.prepared/'weights.npy');trace=K.T@(weights*raw)
    moment=-np.load(checkpoint)['y']
    target=np.zeros(op.free.shape[1]);target[0]=1
    original_ideal_defect=float(np.max(abs(op.free.T@moment-target)))
    sys.path.insert(0,str(ROOT/'.venv-interacting-libs'))
    from research.nvidia_followup_20260915.sparse_quotient import SparseQuotient
    quotient=SparseQuotient(op)
    original_energy=float(op.rhs@moment)
    moment=quotient.y0+quotient.project(moment)
    defect=float(np.max(abs(op.free.T@trace-target)))
    if defect>1e-9:raise ValueError('Singlet trace failed numerical ideal and normalization equations')
    records=[];mixture=0.;null_max=0.
    for index,(physical,proposed) in enumerate(zip(op.AT(trace),op.AT(moment))):
        ev,U=linalg.eigh(physical,check_finite=False)
        if ev[0]<-1e-10:raise ValueError('Singlet trace Gram is not numerically PSD')
        positive=ev>1e-10;null=U[:,~positive]
        null_defect=float(np.max(abs(proposed@null),initial=0))
        null_max=max(null_max,null_defect)
        if positive.any():
            V=U[:,positive]/np.sqrt(ev[positive])[None,:]
            minimum=float(linalg.eigvalsh(V.T@proposed@V,subset_by_index=(0,0),check_finite=False)[0])
            needed=max(0.,-minimum/(1-minimum)) if minimum<0 else 0.
        else:minimum=0.;needed=0.
        mixture=max(mixture,needed)
        records.append({'block':index,'dimension':len(ev),'trace_nullity':int((~positive).sum()),
            'proposed_null_action_max':null_defect,'minimum_relative_positive_face_eigenvalue':minimum,
            'face_only_mixture_needed':needed})
    U=float(Fraction(op.meta['upper_Ha']));E=float(op.rhs@moment);Et=float(op.rhs@trace)
    record={'kind':'approximate_singlet_trace_dual_diagnostic','source':str(checkpoint),
        'seconds':time.monotonic()-started,'singlet_dimension':dimension,'trace_ideal_defect':defect,
        'original_proposal_ideal_defect':original_ideal_defect,'original_proposal_energy_Ha':original_energy,
        'numerical_affine_projection_used':True,'exact_affine_projection_used':False,
        'approximate_dual_energy_Ha':E,'unverified_family_floor_mHa':1000*(U-E),
        'trace_energy_Ha':Et,'maximum_null_action':null_max,
        'positive_face_only_mixture':mixture,'face_only_energy_after_mixing_Ha':(1-mixture)*E+mixture*Et,
        'affine_and_nullspace_repair_performed':False,'exact_PSD_proof_performed':False,
        'family_obstruction_proved':False,'blocks':records,'lift':details}
    with output.open('x') as stream:json.dump(record,stream,indent=2);stream.write('\n')
    print(json.dumps({k:v for k,v in record.items() if k not in ('blocks','lift')}),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('case',type=Path);p.add_argument('checkpoint',type=Path)
    p.add_argument('output',type=Path);a=p.parse_args();run(a.case.resolve(),a.checkpoint.resolve(),a.output.resolve())
