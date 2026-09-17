"""Exact scalar part of a conditional spectral witness.

The caller must prove overlap, vector norm and recurrence residual bounds.
This gate alone does not verify a Hamiltonian or tensor data.
"""
from fractions import Fraction as F


def rational(x):
    if isinstance(x,bool) or not isinstance(x,(int,str,F)):
        raise TypeError('Use integers, rational strings or Fractions, not floats')
    return F(x)


def gate(b,ell,z,norm_upper,errors,seed_error=F(0),overlap=F(1)):
    """Sufficient lower-energy gate using a conservative geometric envelope."""
    b,ell,z,norm_upper,seed_error,overlap=map(rational,(b,ell,z,norm_upper,seed_error,overlap))
    errors=[rational(x) for x in errors]
    if not (0<z<1 and ell<b):raise ValueError('Require 0<z<1 and ell<b')
    if norm_upper<0 or seed_error<0 or any(x<0 for x in errors):
        raise ValueError('Norm and error bounds must be nonnegative')
    if overlap<=seed_error or not errors:
        raise ValueError('Require positive net overlap and at least one step')
    k=len(errors)
    score=2*z**k*norm_upper+seed_error
    score+=2*sum((z**(i+1)*e for i,e in enumerate(errors)),F(0))/(1-z*z)
    lower=(b+ell)/2-(b-ell)*(z+1/z)/4
    return {'accepted':score<overlap,'score':score,'overlap':overlap,'lower':lower,'k':k}
