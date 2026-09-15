"""Exact original-H recognition for a restricted, interacting path parent family.

Accepted H = c + sum_i g_i B_i^dagger B_i, where g_i>0 and
B_i = q_(i+1) a_i (1-n_(i+1)) - q_i a_(i+1) (1-n_i).
The nonzero fixed-N vector with amplitudes product(q_i) is a common null
vector. This proves E0=c without enumerating an occupation sector.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse
import json
import time
from experiments.marginal_symbolic import (add, adj, canonical, decode, encode, hermitian,
    mono, product, scale, validate_word)


def edge_factor(i, left, right):
    ni=mono(((1,i),(0,i)))
    nj=mono(((1,i+1),(0,i+1)))
    return add(scale(product(mono(((0,i),)),add(mono(()),scale(nj,-1))),right),
               scale(product(mono(((0,i+1),)),add(mono(()),scale(ni,-1))),-left))


def _recognize_ordered(poly, modes, particles):
    start=time.monotonic()
    if type(modes) is not int or modes<2 or type(particles) is not int or not 0<=particles<=modes:
        raise ValueError('Invalid sector')
    if not isinstance(poly,dict):raise ValueError('Polynomial must be dict')
    for w,c in poly.items():
        validate_word(w,modes,4)
        if type(c) is bool or not isinstance(c,(int,F)):raise ValueError('Exact coefficients required')
    h=canonical(poly)
    def reject(reason):return {'accepted':False,'reason':reason,'wall_seconds':time.monotonic()-start}
    def n(i):return ((1,i),(0,i))
    A=[F(h.get(n(i),0)) for i in range(modes)]
    t=[-F(h.get(((1,i),(0,i+1)),0)) for i in range(modes-1)]
    # n_i*n_j is MINUS the canonical word with ascending annihilators.
    u=[F(h.get(((1,i),(1,i+1),(0,i),(0,i+1)),0)) for i in range(modes-1)]
    if any(x<=0 for x in t):return reject('Positive nearest-neighbor hopping magnitudes required')
    alpha=[];beta=[];previous=F(0)
    for i in range(modes-1):
        a=A[i]-previous
        if a<=0:return reject('Nonpositive recovered edge diagonal')
        b=t[i]**2/a
        if u[i]!=a+b:return reject('Density coefficient mismatch')
        alpha.append(a);beta.append(b);previous=b
    if A[-1]!=beta[-1]:return reject('Boundary onsite mismatch')
    q=[F(1)];g=[]
    for i in range(modes-1):
        q.append(q[-1]*alpha[i]/t[i]);g.append(t[i]/(q[i]*q[i+1]))
    c=F(h.get((),0));expected=mono((),c);factors=[];pairs=0
    for i in range(modes-1):
        B=edge_factor(i,q[i],q[i+1]);factors.append(B)
        square=scale(product(canonical(adj(B)),B),g[i]);pairs+=len(B)**2
        for w,v in square.items():expected[w]=expected.get(w,F(0))+v
    expected={w:v for w,v in expected.items() if v}
    if expected!=h:return reject('Original Hamiltonian differs from reconstructed positive parent')
    dp=[F(1)]+[F(0)]*particles;updates=0;max_dp_bits=1
    for x in q:
        for k in range(particles,0,-1):
            dp[k]+=x*x*dp[k-1];updates+=1
            max_dp_bits=max(max_dp_bits,abs(dp[k].numerator).bit_length(),dp[k].denominator.bit_length())
    if dp[particles]<=0:raise AssertionError('Common null state must be nonzero')
    returned=[c,*q,*g,*alpha,*beta,dp[particles]]
    receipt={'accepted':True,'modes':modes,'particles':particles,'lower':str(c),'upper':str(c),'width':'0',
      'q':[str(x) for x in q],'g':[str(x) for x in g],'normalization':str(dp[particles]),
      'edge_count':modes-1,'hamiltonian_terms':len(h),'factor_terms':sum(map(len,factors)),
      'reconstruction_word_pairs':pairs,'dp_updates':updates,'dp_storage_entries':len(dp),
      'max_returned_rational_bits':max(max(abs(x.numerator).bit_length(),x.denominator.bit_length()) for x in returned),
      'max_dp_rational_bits':max_dp_bits,'wall_seconds':time.monotonic()-start,
      'scope':'Exact natural-order connected path parent only. Positive rational weights; original-H CAR identity and nonzero common fixed-N null state. No generic chemistry claim.',
      'supplied_factors_used':False,'occupation_states_enumerated':0}
    return receipt


def recognize(poly, modes, particles):
    """Infer path order and real sign gauge from the original H, then certify."""
    start=time.monotonic()
    if type(modes) is not int or modes<2 or type(particles) is not int or not 0<=particles<=modes:
        raise ValueError('Invalid sector')
    if not isinstance(poly,dict):raise ValueError('Polynomial must be dict')
    for w,c in poly.items():
        validate_word(w,modes,4)
        if type(c) is bool or not isinstance(c,(int,F)):raise ValueError('Exact coefficients required')
    h=canonical(poly)
    def reject(reason):return {'accepted':False,'reason':reason,'wall_seconds':time.monotonic()-start}
    if not hermitian(h):return reject('Hermitian Hamiltonian required')
    neighbors=[set() for _ in range(modes)]
    for w in h:
        if len(w)==2 and w[0][0]==1 and w[1][0]==0 and w[0][1]!=w[1][1]:
            u,v=w[0][1],w[1][1];neighbors[u].add(v);neighbors[v].add(u)
    ends=[i for i,ns in enumerate(neighbors) if len(ns)==1]
    if len(ends)!=2 or any(len(ns) not in (1,2) for ns in neighbors):return reject('Hopping graph is not a connected path')
    order=[];seen=set();current=min(ends);gauges={current:1}
    while current not in seen:
        order.append(current);seen.add(current)
        nxt=neighbors[current]-seen
        if not nxt:break
        other=next(iter(nxt));c=h[((1,current),(0,other))]
        gauges[other]=-gauges[current]*(1 if c>0 else -1);current=other
    if len(order)!=modes:return reject('Disconnected hopping graph')
    inverse={old:new for new,old in enumerate(order)}
    def transform(p,mode_map):
        out={}
        for w,c in p.items():
            mapped=[]
            for create,i in w:
                j,sign=mode_map[i];mapped.append((create,j));c*=sign
            out[tuple(mapped)]=c
        return canonical(out)
    mapped=transform(h,{i:(inverse[i],gauges[i]) for i in range(modes)})
    # The same sign on creation and annihilation is an exact real CAR gauge.
    restored=transform(mapped,{i:(old,gauges[old]) for i,old in enumerate(order)})
    if restored!=h:raise AssertionError('Signed mode permutation failed round trip')
    result=_recognize_ordered(mapped,modes,particles)
    result.update(path_order=order,path_gauges=[gauges[i] for i in order],
                  signed_permutation_round_trip_checked=True,wall_seconds=time.monotonic()-start)
    if result['accepted']:
        result['scope']='Exact connected hopping-path parent after inferred signed mode permutation. Original-H CAR identity and transported common fixed-N null state. No generic chemistry claim.'
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--fixture',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();f=json.loads(a.fixture.read_text())
    r=recognize(decode(f['hamiltonian'],f['modes'],4),f['modes'],f['particles'])
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(r,indent=2)+'\n')
    print(json.dumps(r))

if __name__=='__main__':main()
