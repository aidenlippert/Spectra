"""Exact policy-cover certificate kernel; premises are supplied, not inferred.

The bounds concern whatever response the caller's observations actually enclose.
This module cannot validate physical applicability or manufacture a policy decoder.
"""
from fractions import Fraction as F
from dataclasses import dataclass


def rational(x):
    if type(x) not in (int,F):raise ValueError('exact rational required')
    x=F(x)
    if max(abs(x.numerator).bit_length(),x.denominator.bit_length())>4096:
        raise ValueError('rational budget')
    return x


@dataclass(frozen=True)
class Sample:
    point:tuple
    intervals:tuple


def box_for(path,dim):
    if type(dim) is not int or not 1<=dim<=8:raise ValueError('dimension budget')
    if not isinstance(path,str) or len(path)>64 or set(path)-{'0','1'}:raise ValueError('invalid cell path')
    box=[(F(0),F(1)) for _ in range(dim)]
    for bit in path:
        axis=max(range(dim),key=lambda j:box[j][1]-box[j][0])
        lo,hi=box[axis];mid=(lo+hi)/2
        box[axis]=(lo,mid) if bit=='0' else (mid,hi)
    return tuple(box)


def _cover(paths,dim):
    if not isinstance(paths,(list,tuple)) or not 1<=len(paths)<=512:raise ValueError('cell budget')
    if len(set(paths))!=len(paths):raise ValueError('duplicate cells')
    boxes={p:box_for(p,dim) for p in paths}
    ordered=sorted(paths)
    if any(b.startswith(a) for a,b in zip(ordered,ordered[1:])):raise ValueError('overlapping tree leaves')
    if sum((F(1,2**len(p)) for p in paths),F(0))!=1:raise ValueError('incomplete policy cover')
    return boxes


def _samples(samples,dim,L):
    if not isinstance(samples,(list,tuple)) or not 1<=len(samples)<=512:raise ValueError('observation budget')
    if not isinstance(L,(list,tuple)) or not 1<=len(L)<=16:raise ValueError('quantity budget')
    L=tuple(rational(v) for v in L)
    if min(L)<0:raise ValueError('negative modulus')
    normalized=[]
    for s in samples:
        if not isinstance(s,Sample) or len(s.point)!=dim or len(s.intervals)!=len(L):raise ValueError('sample shape')
        point=tuple(rational(v) for v in s.point)
        if any(not 0<=v<=1 for v in point):raise ValueError('point outside domain')
        intervals=[]
        for pair in s.intervals:
            if not isinstance(pair,(list,tuple)) or len(pair)!=2:raise ValueError('interval shape')
            lo,hi=map(rational,pair)
            if lo>hi:raise ValueError('reversed interval')
            intervals.append((lo,hi))
        normalized.append(Sample(point,tuple(intervals)))
    return normalized,L


def _bounds(box,samples,L):
    result=[]
    for j,lipschitz in enumerate(L):
        lows=[];highs=[]
        for s in samples:
            distance=max(max(abs(x-lo),abs(x-hi)) for x,(lo,hi) in zip(s.point,box))
            lows.append(s.intervals[j][0]-lipschitz*distance)
            highs.append(s.intervals[j][1]+lipschitz*distance)
        result.append((max(lows),min(highs)))
    return tuple(result)


def certify(paths,samples,L,dim,epsilon):
    """Replay a complete binary partition and exact Lipschitz implications."""
    epsilon=rational(epsilon)
    if epsilon<0:raise ValueError('negative requested gap')
    boxes=_cover(paths,dim);samples,L=_samples(samples,dim,L)
    points=[_bounds(tuple((x,x) for x in s.point),samples,L) for s in samples]
    if any(lo>hi for bs in points for lo,hi in bs):
        return dict(status='inconsistent_assumptions',reason='observations and modulus conflict')
    cells={p:_bounds(b,samples,L) for p,b in boxes.items()}
    active=[p for p,bs in cells.items() if all(lo<=0 for lo,_ in bs[1:])]
    if not active:return dict(status='no_feasible_policy_in_covered_domain',cells=cells)
    lower=min(cells[p][0][0] for p in active)
    feasible=[i for i,bs in enumerate(points) if all(hi<=0 for _,hi in bs[1:])]
    if not feasible:return dict(status='unresolved',lower=lower,upper=None,active=active,cells=cells)
    best=min(feasible,key=lambda i:points[i][0][1]);upper=points[best][0][1]
    if upper<lower:raise AssertionError('invalid certificate bracket')
    return dict(status='conditional_certificate' if upper-lower<=epsilon else 'unresolved',
        lower=lower,upper=upper,gap=upper-lower,policy=samples[best].point,
        policy_bounds=points[best],active=active,cells=cells)


def construct(oracle,L,dim,epsilon,max_observations=96):
    """Bounded adaptive cover refinement; no statistical oracle is hidden here.

    oracle(point,requested_radius) must return justified response intervals.
    The request is advisory: wider returned intervals retain their full width.
    Each call is charged. The returned policy is the declared coordinate word.
    """
    if type(max_observations) is not int or not 1<=max_observations<=256:raise ValueError('observation budget')
    box_for('',dim);L=tuple(L);epsilon=rational(epsilon)
    paths=[''];samples=[];chosen='';trace=[]
    for iteration in range(max_observations):
        box=box_for(chosen,dim);point=tuple((lo+hi)/2 for lo,hi in box)
        requested=F(1,2**(len(chosen)+3))
        intervals=oracle(point,requested)
        samples.append(Sample(point,tuple(intervals)))
        receipt=certify(paths,samples,L,dim,epsilon)
        trace.append(dict(cell=chosen,point=point,requested_radius=requested,status=receipt['status']))
        if receipt['status']!='unresolved':break
        if iteration+1==max_observations:break
        chosen=min(receipt['active'],key=lambda p:(receipt['cells'][p][0][0],len(p),p))
        if len(chosen)>=63:break
        paths.remove(chosen);paths.extend([chosen+'0',chosen+'1'])
        # Observe the less costly child first; both remain in the complete cover.
        chosen=min((chosen+'0',chosen+'1'),key=lambda p:(_bounds(box_for(p,dim),samples,L)[0][0],p))
    return dict(receipt=receipt,paths=paths,samples=samples,trace=trace,
                observation_calls=len(samples),leaf_cells=len(paths),
                scope='conditional interval arithmetic; no physical validation')
