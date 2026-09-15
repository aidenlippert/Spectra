"""Exact SU(2) Haar projector on even CAR polynomials of degree at most six."""
from fractions import Fraction as F
from experiments.marginal_symbolic import canonical,mono,add,scale
from research.certificate_scaling.spin_basis import spin_weight2


def ladder(poly,raising):
    out={}
    for w,v in poly.items():
        for k,(c,i) in enumerate(w):
            if not ((c==1 and i%2==(1 if raising else 0)) or (c==0 and i%2==(0 if raising else 1))):continue
            word=w[:k]+((c,i^1),)+w[k+1:];factor=v if c else -v
            for u,s in canonical(mono(word)).items():out[u]=out.get(u,F(0))+factor*s
    return {w:v for w,v in out.items() if v}


def casimir(poly):
    z2={w:v*F(spin_weight2(w)**2,4) for w,v in poly.items() if spin_weight2(w)}
    return add(z2,scale(add(ladder(ladder(poly,False),True),ladder(ladder(poly,True),False)),F(1,2)))


def twirl(poly):
    if any(len(w)>6 or len(w)%2 for w in poly):raise ValueError('Even degree at most six required')
    if any(type(c) is bool or not isinstance(c,(int,F)) for c in poly.values()):raise ValueError('Exact real coefficients required')
    out=canonical(poly);jmax=max((len(w)//2 for w in out),default=0)
    for j in range(1,jmax+1):out=add(out,scale(casimir(out),F(-1,j*(j+1))))
    if casimir(out):raise AssertionError('Twirl failed exact Casimir kernel check')
    return out
