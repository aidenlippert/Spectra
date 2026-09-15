"""CAR coupling bound for P={N_high<=1} to R={N_high>=2}."""
from fractions import Fraction as F
import time
from math import isqrt
from experiments.marginal_symbolic import decode, canonical, add, word_product
from experiments.marginal_hunt_car import adj
from experiments.marginal_transfer_verify import apply_word
from experiments.marginal_general_schur import norm_bound

def bounded_product(p,q,deadline,stats):
    out={}
    for left,a in p.items():
        for right,b in q.items():
            stats['word_products']+=1
            if time.monotonic()>deadline: raise ValueError('Operator construction time cap exceeded')
            for w,sgn in word_product(left,right): out[w]=out.get(w,F(0))+a*b*sgn
            if len(out)>100000: raise ValueError('Canonical operator word cap exceeded')
    return {w:c for w,c in out.items() if c}

def sqrt_up(v,den=10**8):
    if v<0: raise ValueError('Negative squared norm')
    q=isqrt(v.numerator*den*den//v.denominator)
    if F(q*q,den*den)<v:q+=1
    return F(q,den)

def local_blocks(data):
    m=data['modes']; r=m-4
    if m not in (12,16) or data['particles']!=m//2: raise ValueError('H6/H8 required')
    inputs=[0,1,2,4,8]; outputs=[q for q in range(16) if q.bit_count()>=2]
    polys=[[{} for _ in inputs] for _ in outputs]; h=decode(data['hamiltonian'],m,4)
    for word,coefficient in h.items():
        rest=tuple((c,i) for c,i in word if i<r); local=tuple((c,i-r) for c,i in word if i>=r)
        crossings=sum(i>=r and j<r for a,(_,i) in enumerate(word) for _,j in word[a+1:])
        for col,source in enumerate(inputs):
            value=apply_word(local,source)
            if value and value[0] in outputs:
                row=outputs.index(value[0]); parity=crossings+len(local)*(data['particles']-source.bit_count())
                polys[row][col][rest]=polys[row][col].get(rest,F(0))+coefficient*value[1]*(-1 if parity%2 else 1)
    return [[canonical(p) for p in row] for row in polys],inputs,outputs

def prove(data):
    start=time.monotonic(); blocks,inputs,outputs=local_blocks(data); deadline=start+60; stats={'word_products':0}; norms=[]; rows=[]
    for row in blocks:
        nr=[]
        for p in row:
            if not p: nr.append(F(0)); continue
            degrees={len(w)%2 for w in p}
            if len(degrees)!=1: raise ValueError('Mixed fermion parity')
            sq=bounded_product(adj(p),p,deadline,stats)
            if next(iter(degrees)): sq=add(sq,bounded_product(p,adj(p),deadline,stats))
            mu=norm_bound(sq,'hermitian_pairs'); nr.append(sqrt_up(mu)); rows.append({'terms':len(p),'squared_norm_bound':str(mu)})
        norms.append(nr)
    nu=max(map(sum,norms))*max(map(sum,zip(*norms)))
    return {'coupling_norm_squared_Ha2':str(nu),'block_norm_bounds':[[str(x) for x in row] for row in norms],
      'local_input_occupations':inputs,'local_output_occupations':outputs,'local_occupation_labels':len(inputs)+len(outputs),
      'local_modes':4,'remaining_modes_not_enumerated':data['modes']-4,'blocks':rows,
      'many_body_states_enumerated':0,'many_body_matrix_entries':0,'word_products':stats['word_products'],
      'wall_seconds':time.monotonic()-start}
