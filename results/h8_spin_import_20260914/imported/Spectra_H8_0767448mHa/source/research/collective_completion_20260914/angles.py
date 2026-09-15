"""Audited filter mathematics and bounded molecular transition probes.

These diagnostics do not accept an energy bound. They explicitly count the
few determinant labels queried; no complete determinant basis is generated.
"""
from fractions import Fraction as F
from itertools import combinations
import json
from math import comb,sqrt,acosh,ceil
from pathlib import Path
from experiments.marginal_symbolic import decode
from research.reconstruction_compression_20260914.inputs import dump,sha

ROOT=Path(__file__).resolve().parents[2];OUT=ROOT/'results/collective_completion_20260914'

def act(word,state):
    sign=1
    for c,i in reversed(word):
        bit=1<<i
        if bool(state&bit)==bool(c):return None,0
        if (state&(bit-1)).bit_count()%2:sign=-sign
        state^=bit
    return state,sign

def neighbors(h,state):
    out={}
    for w,c in h.items():
        x,s=act(w,state)
        if s:out[x]=out.get(x,F(0))+s*c
    return {x:c for x,c in out.items() if c and x!=state}

def amplitude(cert,state):
    E={0:1}
    for i,edges in enumerate(cert['tensors']):
        out={}
        for a,s,b,v in edges:
            if s==((state>>i)&1) and a in E:out[b]=out.get(b,0)+E[a]*v
        E={b:v for b,v in out.items() if v}
    return E.get(0,0)

def cheb(kind,k,x):
    if k<0:return F(0)
    a,b=F(1),x if kind=='T' else 2*x
    if k==0:return a
    for j in range(1,k):a,b=b,2*x*b-a
    return b

def filter_audit():
    # The exact coefficients of a defect R_j, relative to T_k(x), are
    # T_j(x_b) U_{k-j}(x)/T_k(x). Check their squared universal bound using
    # rational arithmetic at many points; the proof is in MATHEMATICS.md.
    checked=0;largest=0.
    for xb in (F(1001,1000),F(11,10),F(3,2),F(2)):
        for x in (xb,xb+F(1,100),xb+1,xb+10):
            for k in (1,2,3,8,24,48):
                for j in range(1,k+1):
                    w=cheb('T',j,xb)*cheb('U',k-j,x)/cheb('T',k,x)
                    assert w*w<=xb*xb/(xb*xb-1)
                    checked+=1;largest=max(largest,float(w))
    return {'exact_scalar_kernel_checks':checked,'status':'derived sufficient theorem, no molecular filter witness',
            'supplied_toy_archives_available_locally':False,'molecular_certificate_accepted':False,
            'required_final_norm':'Unnormalized Hilbert-Schmidt norm on the whole target sector',
            'missing_construction':'Compact filter sequence plus deterministic final trace and all operator-norm defect bounds'}

def molecular(case):
    cfg=json.loads((OUT/'frozen_inputs.json').read_text())[case];data=json.loads((ROOT/cfg['fixture']).read_text());state=json.loads((ROOT/cfg['state']).read_text())
    h=decode(data['hamiltonian'],data['modes'],4);m,n=data['modes'],data['particles'];root=(1<<n)-1
    edges=neighbors(h,root);selected=sorted(edges,key=lambda x:(-abs(edges[x]),x))[:48];cycle=None
    # Only edges between the chosen neighbors are queried. Stop on the first
    # exact frustrated triangle through the reference occupation pattern.
    labels={root,*selected};queries=1
    for a in selected:
        around=neighbors(h,a);queries+=1
        for b in selected:
            if b<=a or b not in around:continue
            if edges[a]*around[b]*edges[b]>0:
                cycle={'states':[root,a,b],'edge_Ha':[str(edges[a]),str(around[b]),str(edges[b])],
                       'cycle_product_positive':True,'all_negative_diagonal_sign_gauge_impossible_for_this_cycle':True};break
        if cycle:break
    amps=[amplitude(state,x) for x in sorted(labels)]
    # A spectral enclosure from the coefficient norm is deliberately crude,
    # but verified. Its scale exposes the filter degree/precision demand.
    constant=h.get((),F(0));radius=sum(abs(c) for w,c in h.items() if w);upper_spectrum=constant+radius
    U=F(cfg['upper_Ha']);b=U-F(1,625);a=b+F(1,1250);xb=(a+upper_spectrum-2*b)/(upper_spectrum-a)
    D=comb(m,n);density_threshold=F(1,D)
    # This degree suppresses the interval [a,c] uniformly. It is a diagnostic,
    # not a lower proof: no assumption is made that the spectrum lies above a.
    degree=ceil(acosh(sqrt(D))/acosh(float(xb)))
    return {'case':case,'fixture_sha256':sha(ROOT/cfg['fixture']),'state_sha256':sha(ROOT/cfg['state']),
            'complete_sector_dimension':D,'selected_determinant_labels':len(labels),'Hamiltonian_column_queries':queries,
            'full_determinant_basis_generated':False,'signed_cycle':cycle,
            'selected_MPS_amplitude_signs':{'positive':sum(v>0 for v in amps),'negative':sum(v<0 for v in amps),'exact_zero':sum(v==0 for v in amps)},
            'longest_nonconstant_operator_span_spin_modes':max(max(i for c,i in w)-min(i for c,i in w)+1 for w in h if w),
            'spectral_filter_diagnostic':{'target_b_Ha':str(b),'a_Ha':str(a),'certified_spectral_upper_Ha':str(upper_spectrum),
                'normalized_squared_trace_must_be_below':str(density_threshold),'rough_uniform_suppression_degree':degree,
                'defect_multiplier_approx':float(xb)/sqrt(float(xb*xb-1)),'accepts_a_molecular_bound':False},
            'scope':'A cycle excludes only diagonal sign gauges in this occupation basis. It does not exclude signed routing, another basis, or another compact ansatz.'}

if __name__=='__main__':
    receipt={'filter':filter_audit(),'molecules':{c:molecular(c) for c in ('h6','h8')}}
    dump(OUT/'alternative_angles.json',receipt);print(json.dumps(receipt,indent=2))
