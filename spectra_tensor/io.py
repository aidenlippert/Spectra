"""Portable Gaussian-integer tensor witnesses."""
from pathlib import Path
from fractions import Fraction as F
import json
import numpy as np
from . import tensor as t,exact

def export_program(x,bits=30):
    if type(bits) is not int or not 1<=bits<=40:raise ValueError('rounding bits')
    x.validate();scale=1<<bits;cores=[]
    for a in x.cores:
        real=np.rint(a.real*scale);imag=np.rint(a.imag*scale)
        if np.any(abs(real)>=1<<48) or np.any(abs(imag)>=1<<48):raise ValueError('coefficient envelope')
        real=real.astype(np.int64);imag=imag.astype(np.int64)
        cores.append([[int(l),int(p),int(r),int(real[l,p,r]),int(imag[l,p,r])] for l,p,r in zip(*np.nonzero((real!=0)|(imag!=0)))])
    return dict(kind='charge_mps_gaussian_integer_v1',local_denominator=scale,charges=[[list(q) for q in cut] for cut in x.charges],cores=cores)

def import_program(program):
    n=len(program['cores']);exact.validate_program(program,n);qs=[[tuple(q) for q in cut] for cut in program['charges']];cores=[];scale=program['local_denominator']
    for i,entries in enumerate(program['cores']):
        a=np.zeros((len(qs[i]),4,len(qs[i+1])),complex)
        for l,p,r,re,im in entries:a[l,p,r]=(re+1j*im)/scale
        cores.append(a)
    return t.MPS(cores,qs).validate()

def request(spec,zs,target='1/1000'):
    return dict(kind='tensor_hubbard_response_request_v1',model=spec,source=exact.source(spec),frequencies=[dict(omega=str(F(float(z.real))),eta=str(F(float(z.imag)))) for z in zs],target_radius=target)

def candidate(spec,zs,xs,bits=30):
    if len(zs)!=len(xs):raise ValueError('frequency count')
    return dict(kind='tensor_hubbard_response_v1',model_hash=exact.model_hash(spec),queries=[dict(omega=str(F(float(z.real))),eta=str(F(float(z.imag))),program=export_program(x,bits)) for z,x in zip(zs,xs)])

def dump(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('x') as f:json.dump(value,f,separators=(',',':'))
