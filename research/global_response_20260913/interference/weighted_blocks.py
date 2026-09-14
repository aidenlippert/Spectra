"""Weighted block-Schur domination for the local fermionic coupling.

The existing coupling proof uses max row sum times max column sum.  This
retains each occupation-conditioned anticommutator block but permits a
positive diagonal metric on the source and target occupation labels.  For
nonnegative block norm bounds n_ij, the induced weighted row/column bounds
give ||B||^2 <= R(w) C(w), an exact sufficient rule for every positive
weight vector.  We only use powers-of-two weights for replayable diagnostics.
"""
from fractions import Fraction as F
from itertools import product
from research.composable_response_20260913 import coupling

def weighted_bound(norms, exponents):
    if len(exponents) != len(norms)+len(norms[0]): raise ValueError('weight dimension mismatch')
    if any(type(e) is not int or e < -32 or e > 32 for e in exponents): raise ValueError('invalid weight exponent')
    r=len(norms); c=len(norms[0]); wr=[F(2)**e for e in exponents[:r]]; wc=[F(2)**e for e in exponents[r:]]
    rows=[sum(norms[i][j]*wc[j]/wr[i] for j in range(c)) for i in range(r)]
    cols=[sum(norms[i][j]*wr[i]/wc[j] for i in range(r)) for j in range(c)]
    return max(rows)*max(cols), max(rows), max(cols)

def optimize(data, radius=1):
    if type(radius) is not int or radius < 0 or radius > 1:
        raise ValueError('bounded weight search requires radius 0 or 1')
    p=coupling.prove(data); n=[[F(x) for x in row] for row in p['block_norm_bounds']]
    base=max(map(sum,n))*max(map(sum,zip(*n)))
    best=(base,(0,)*(len(n)+len(n[0])))
    for ex in product(range(-radius,radius+1), repeat=len(n)+len(n[0])):
        v,_,_=weighted_bound(n,ex)
        if v<best[0]: best=(v,ex)
    ex=best[1]
    weights=[str(F(2)**e) for e in ex]
    return {'unweighted_squared_norm':str(base),'weighted_squared_norm':str(best[0]),'exponents':list(ex),'weights':weights,
            'improvement_ratio':float(best[0]/base) if base else 1.0,'coupling':p,
            'rule':'weighted block row/column Schur bound over anticommutator blocks'}
