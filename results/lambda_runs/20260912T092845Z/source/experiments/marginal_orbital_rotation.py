"""Exact spatial-orbital rotation probe for the H6 marginal geometry."""
from fractions import Fraction as F
from itertools import combinations
import json, math
from pathlib import Path

from experiments.marginal_hunt_car import add, adj, mono, mul
from experiments.marginal_symbolic import decode, encode
from experiments.marginal_determinant_tree import DeterminantOracle
from experiments.marginal_spin_reduction import SpinZeroOracle
from experiments.marginal_spin_constructor import select_reference, spin_states

def rotation(modes, left=2, right=3, c=F(3,5), s=F(4,5)):
    if any(type(x) is not int for x in (modes,left,right)) or modes % 2 or not 0 <= left < modes//2 or not 0 <= right < modes//2 or left == right:
        raise ValueError('Invalid spatial rotation')
    if not isinstance(c,F) or not isinstance(s,F): raise ValueError('Rotation entries must be exact rationals')
    if c*c+s*s != 1: raise ValueError('Rotation is not orthogonal')
    u = [[F(int(i==j)) for j in range(modes//2)] for i in range(modes//2)]
    u[left][left]=u[right][right]=c; u[left][right]=s; u[right][left]=-s
    return u

def rotate(poly, u):
    if not isinstance(u,list) or not u or any(len(row)!=len(u) or any(not isinstance(x,F) for x in row) for row in u):
        raise ValueError('Rotation matrix must be rational and square')
    if any(sum(u[i][k]*u[j][k] for k in range(len(u))) != (1 if i==j else 0) for i in range(len(u)) for j in range(len(u))):
        raise ValueError('Rotation matrix must be orthogonal')
    n=2*len(u); out={}
    for word, coeff in poly.items():
        terms={():coeff}
        for cr, mode in word:
            orb, spin=divmod(mode,2)
            image={}
            for w,a in terms.items():
                for q in range(len(u)):
                    if not u[orb][q]: continue
                    image[w+((cr,2*q+spin),)] = image.get(w+((cr,2*q+spin),),F(0))+a*u[orb][q]
            terms=image
        # product normal-orders and combines CAR signs
        expanded={}
        for w,a in terms.items():
            for z,b in mul(mono(w), mono(())).items(): expanded[z]=expanded.get(z,F(0))+a*b
        for w,a in expanded.items(): out[w]=out.get(w,F(0))+a
    return {w:a for w,a in out.items() if a}

def car_check(u):
    n=2*len(u)
    # Direct coefficient check for {a'_p,a'^dag_q}=delta_pq.
    for p in range(n):
        op,sp=divmod(p,2)
        for q in range(n):
            oq,sq=divmod(q,2)
            value=(sum(u[op][k]*u[oq][k] for k in range(len(u))) if sp==sq else F(0))
            if value != (1 if p==q else 0): return False
    return True

def sector_spectrum(data):
    oracle=SpinZeroOracle(data)
    states=sorted(spin_states(oracle)); idx={s:i for i,s in enumerate(states)}
    import numpy as np
    a=np.zeros((len(states),len(states)))
    for j,s in enumerate(states):
        for t,v in oracle.action(s).items(): a[idx[t],j]=float(v)
    return np.linalg.eigvalsh(a)

def run(source=None, output=None):
    root=Path(__file__).resolve().parents[1]
    source=Path(source or root/'results/marginal_h6/direct_spin/symmetric_hamiltonian.json')
    d=json.loads(source.read_text()); h=decode(d['hamiltonian'],d['modes'],4); u=rotation(d['modes'])
    if not car_check(u): raise AssertionError('CAR check failed')
    hr=rotate(h,u); rotated=dict(d,hamiltonian=encode(hr))
    import numpy as np
    e0=sector_spectrum(d); e1=sector_spectrum(rotated)
    def q_data(x):
        o=SpinZeroOracle(x); vals=[]
        px,_,_=select_reference(o,32)
        remaining=set(spin_states(o))-o.retained(px)
        physical=[]; comparison=[]; dims=[]
        while remaining:
            seed=min(remaining); group={seed}; remaining.remove(seed); todo=[seed]
            while todo:
                z=todo.pop()
                for y in o.action(z):
                    if y in remaining: remaining.remove(y); group.add(y); todo.append(y)
            mat=np.array([[float(o.action(a).get(b,0)) for b in sorted(group)] for a in sorted(group)])
            physical.append(float(np.linalg.eigvalsh(mat)[0])); dims.append(len(group))
            comparison.append(float(np.linalg.eigvalsh(np.diag(np.diag(mat))-np.abs(mat-np.diag(np.diag(mat))))[0]))
        return {'p_states':px,'block_dimensions':dims,'physical_floor':min(physical),'comparison_floor':min(comparison)}
    q0=q_data(d); q1=q_data(rotated)
    result={'rotation':{'left':2,'right':3,'cos':'3/5','sin':'4/5'},'car_valid':True,
            'sz0_dimension':len(e0),'spectrum_max_abs_error':float(np.max(np.abs(e0-e1))),
            'original_q':q0,'rotated_q':q1,
            'comparison_floor_delta':q1['comparison_floor']-q0['comparison_floor'],
            'physical_floor_delta':q1['physical_floor']-q0['physical_floor'],
            'inverse_exact':rotate(hr,rotation(d['modes'],c=F(3,5),s=-F(4,5)))==h,
            'scope':'Finite H6 diagnostic; no certificate or scalability claim.'}
    if output: Path(output).write_text(json.dumps(result,indent=2)+'\n')
    return result

if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser(); ap.add_argument('--source'); ap.add_argument('--output'); args=ap.parse_args()
    print(json.dumps(run(args.source,args.output),indent=2))
