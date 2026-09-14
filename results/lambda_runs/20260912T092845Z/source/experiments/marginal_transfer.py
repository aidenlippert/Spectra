"""Connected six-mode transfer test for coefficient-space certificates."""
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import json, argparse
import numpy as np

from experiments.marginal_symbolic import add, mono
from experiments.marginal_coefficient import dictionaries, solve_coefficients, export

def connected_hamiltonian(delta, t=F(1,5)):
    # Intragroup repulsion; matched hopping plus a weak edge joining the two
    # components left by the asymmetric control graph.
    h={}
    for group in ((0,1,2),(3,4,5)):
        for i,j in combinations(group,2): h=add(h,mono(((1,i),(0,i),(1,j),(0,j)),1))
    edges=((0,3,F(1)),(1,4,F(7,10)),(2,5,F(13,10)),(0,4,F(2,5)),(1,2,F(delta)))
    for i,j,w in edges:
        h=add(h,mono(((1,i),(0,j)),-t*w),mono(((1,j),(0,i)),-t*w))
    return h

def apply_word(word, state, modes):
    amp=1; bits=state
    for creation, i in reversed(word):
        sign=(-1)**((bits & ((1<<i)-1)).bit_count())
        if creation:
            if bits>>i&1:return None
            bits|=1<<i
        else:
            if not(bits>>i&1):return None
            bits&=~(1<<i)
        amp*=sign
    return bits,amp

def upper_ed(h, modes, particles):
    states=[s for s in range(1<<modes) if s.bit_count()==particles]; pos={s:i for i,s in enumerate(states)}
    mat=np.zeros((len(states),len(states)))
    for w,c in h.items():
        for s in states:
            z=apply_word(w,s,modes)
            if z: mat[pos[z[0]],pos[s]]+=float(c)*z[1]
    vals,vecs=np.linalg.eigh(mat); v=vecs[:,0]; scale=10**8
    amps=[int(round(x*scale)) for x in v]
    norm=sum(x*x for x in amps); energy=0
    for w,c in h.items():
        for s,a in zip(states,amps):
            z=apply_word(w,s,modes)
            if z: energy+=F(c)*a*amps[pos[z[0]]]*z[1]
    return {"upper":str(F(energy,norm)),"upper_float":float(F(energy,norm)),"norm":str(norm),"amplitudes":amps,"dimension":len(states)}

def run(delta, family="mixed"):
    h=connected_hamiltonian(F(delta)); blocks=dictionaries(6,family)
    sol=solve_coefficients(h,6,3,blocks); cert,receipt=export(h,6,3,blocks,sol)
    upper=upper_ed(h,6,3); cert["independent_upper"]=upper
    receipt.update({"delta":str(F(delta)),"family":family,"upper":upper["upper"],"width":str(F(upper["upper"])-F(receipt["lower"])),"width_float":float(F(upper["upper"])-F(receipt["lower"]))})
    return cert,receipt

if __name__=="__main__":
    p=argparse.ArgumentParser(); p.add_argument("--delta",default="1/10"); p.add_argument("--family",default="mixed"); a=p.parse_args()
    root=Path(__file__).resolve().parents[1]; out=root/"results/marginal_transfer"; out.mkdir(exist_ok=True)
    c,r=run(a.delta,a.family); name=f"m6_delta{str(F(a.delta)).replace('/','_')}_{a.family}"
    (out/(name+".json")).write_text(json.dumps(c)+"\n"); (out/(name+"_receipt.json")).write_text(json.dumps(r,indent=2)+"\n"); print(json.dumps(r))
