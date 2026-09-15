"""Shared spatial contractions generate coupled spin operators on each fixture."""
from fractions import Fraction as F
import time
import numpy as np

from research.molecular_collective_20260913.core import factor_operators
from research.spin_completion_20260913.discovery import DIRECTION_DENOMINATOR


def generate(model, powers):
    start = time.monotonic()
    if tuple(powers) not in ((0.,), (.5,), (1.,), (0., .5), (0., 1.), (.5, 1.)):
        raise ValueError('Unregistered contraction powers')
    spatial = model.p['spatial']; factors = factor_operators(model.p, model.tail)
    weights = np.array([float(w) for w, _ in factors])
    matrices = np.array([[[float(q.get(((1,2*i),(0,2*j)), 0))
        for j in range(spatial)] for i in range(spatial)] for _, q in factors])
    lookup = {tuple(ref): (gid, j) for gid, frame in enumerate(model.frames)
        for j, ref in enumerate(frame['generators'])}
    candidates = []
    for mode in range(model.p['modes']):
        gid, j = lookup[(-1,-1,-1,mode)]
        v = np.zeros(len(model.frames[gid]['generators'])); v[j] = 1
        candidates.append((gid,v,{'linear_mode': mode}))
    for alpha in powers:
        for i in range(spatial):
            for s in range(2):
                for t in range(2):
                    for u in range(2):
                        parts = {}
                        for k in range(len(factors)):
                            for p in range(spatial):
                                c = weights[k]**alpha*matrices[k,i,p]
                                if c:
                                    gid, j = lookup[(k,s,t,2*p+u)]
                                    v = parts.setdefault(gid,np.zeros(len(model.frames[gid]['generators'])))
                                    v[j] += c
                        for gid,v in sorted(parts.items()):
                            candidates.append((gid,v,{'alpha':alpha,'orbital':i,'spins':[s,t,u]}))
    span = []; diagnostics = []; bases = [[] for _ in model.frames]
    for gid, v, label in candidates:
        v = v/np.linalg.norm(v)
        z = np.rint(v*DIRECTION_DENOMINATOR).astype(np.int64)
        column = model.frame_coefficients[gid]@(z.astype(float)/DIRECTION_DENOMINATOR)
        norm = np.linalg.norm(column); residual = column.copy()
        if bases[gid]:
            Q = np.column_stack(bases[gid])
            for _ in range(2): residual -= Q@(Q.T@residual)
        relative = float(np.linalg.norm(residual)/norm) if norm else 0.
        if norm == 0 or relative <= 1e-9:
            diagnostics.append({'group':gid,'template':label,'reason':'numerically redundant operator column',
                'relative_independence':relative}); continue
        bases[gid].append(residual/np.linalg.norm(residual))
        span.append({'group':gid,'vector':list(map(int,z)),'template':label})
    return span, {'seconds':time.monotonic()-start,'powers':list(powers),
        'raw_symmetry_components':len(candidates),'accepted_directions':len(span),
        'dimensions':[len(x) for x in bases],'removed':diagnostics,
        'scope':'A recorded numerical template subspace, not an exact completeness or energy claim.'}
