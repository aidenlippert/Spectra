"""Finite two-qubit source with hidden GF(2) calibration maps."""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction
import numpy as np
AXES={0:np.array([[0,1],[1,0]],complex),1:np.array([[1,0],[0,-1]],complex)}
def density(a,b,sign,visibility=Fraction(1,2)):
    if type(a) is not int or type(b) is not int or type(sign) is not int or a not in AXES or b not in AXES or sign not in (-1,1): raise ValueError("invalid source")
    if not isinstance(visibility,Fraction) or not 0<=visibility<=1: raise ValueError("invalid visibility")
    return (np.eye(4,dtype=complex)+sign*float(visibility)*np.kron(AXES[a],AXES[b]))/4
def born_probability(a,b,sign,action,visibility=Fraction(1,2)):
    if type(a) is not int or type(b) is not int or type(sign) is not int or a not in AXES or b not in AXES or sign not in (-1,1) or type(action) is not tuple or len(action)!=2 or any(type(x) is not int or x not in (0,1) for x in action): raise ValueError("invalid Born input")
    if not isinstance(visibility,Fraction) or not 0<=visibility<=1: raise ValueError("invalid visibility")
    i,j=action
    return Fraction(1,2)*(1+sign*visibility*(i==a)*(j==b))
def _sample(p,n,rng):
    if type(n) is not int or n<0: raise ValueError("shots must be nonnegative integer")
    if p.denominator > np.iinfo(np.int64).max: raise ValueError("probability denominator exceeds exact sampler range")
    return rng.integers(0,p.denominator,size=n)<p.numerator
@dataclass
class World:
    mask_a: tuple[int,...]
    mask_b: tuple[int,...]
    visibility: Fraction=Fraction(1,2)
    readout_error: Fraction=Fraction(1,10)
    def __post_init__(self):
        if type(self.mask_a) is not tuple or type(self.mask_b) is not tuple or not self.mask_a or len(self.mask_a)!=len(self.mask_b) or len(self.mask_a)>16: raise ValueError("masks must have equal length 1..16")
        if any(type(x) is not int or x not in (0,1) for x in self.mask_a+self.mask_b): raise ValueError("masks must be strict integer bits")
        if not isinstance(self.visibility,Fraction) or not 0<=self.visibility<=1: raise ValueError("invalid visibility")
        if not isinstance(self.readout_error,Fraction) or not 0<=self.readout_error<=Fraction(1,2): raise ValueError("invalid readout error")
        self._mask_a=tuple(self.mask_a); self._mask_b=tuple(self.mask_b)
    @staticmethod
    def _parity(mask,context):
        if type(context) is not tuple or len(context)!=len(mask) or any(type(x) is not int or x not in (0,1) for x in context): raise ValueError("invalid context")
        return sum(a*b for a,b in zip(mask,context))%2
    def calibration(self,side,context,shots,rng):
        if side not in ("a","b"): raise ValueError("invalid side")
        mask=self._mask_a if side=="a" else self._mask_b
        bit=self._parity(mask,context); out=_sample(Fraction(bit),shots,rng); flips=_sample(self.readout_error,shots,rng)
        return tuple(int(x) for x in np.logical_xor(out,flips))
    def target(self,context_a,context_b,sign,action,shots,rng):
        a=self._parity(self._mask_a,context_a); b=self._parity(self._mask_b,context_b)
        plus=_sample(born_probability(a,b,sign,action,self.visibility),shots,rng)
        return tuple(1 if x else -1 for x in plus)
SignedHiddenMap=World
