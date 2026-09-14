"""Compact fermionic frame bound.

For a Hermitian CAR polynomial H=sum c_w W_w, split scalar terms and use
||W_w||<=1.  The resulting lower bound is exact-rational and never builds a
many-body basis.  It is deliberately a diagnostic: coefficientwise domination
cannot exploit cancellation between interference paths.
"""
from fractions import Fraction
import json, hashlib

def _f(x): return Fraction(x)

def digest(x):
    return hashlib.sha256(json.dumps(x, sort_keys=True, separators=(',',':')).encode()).hexdigest()

def frame_lower(hamiltonian):
    """Return c0 - sum |c| for non-scalar CAR words, with exact arithmetic."""
    scalar = Fraction(0); radius = Fraction(0); terms = 0
    for item in hamiltonian:
        word = item['word']; c = _f(item['coefficient'])
        if len(word) == 0:
            scalar += c
        else:
            radius += abs(c); terms += 1
    return {'lower': scalar-radius, 'scalar': scalar, 'radius': radius,
            'non_scalar_terms': terms, 'rule':'scalar minus coefficientwise CAR frame radius'}

def certify(hamiltonian, claimed):
    r = frame_lower(hamiltonian)
    b = _f(claimed)
    if r['lower'] < b:
        raise ValueError('frame bound does not reach claimed target')
    return {**{k:str(v) if isinstance(v,Fraction) else v for k,v in r.items()},
            'claimed':str(b), 'accepted':True}
