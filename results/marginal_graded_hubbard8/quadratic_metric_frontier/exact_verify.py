"""Standalone standard-library verification of exact_certificate.json."""
import json, hashlib
from fractions import Fraction as F
from itertools import combinations, product
from math import comb, prod
from pathlib import Path

def require(condition, message='Exact dual verification gate failed'):
    if not condition: raise ValueError(message)

c = json.loads((Path(__file__).parent / "exact_certificate.json").read_text())
require(c["sites"] == 8 and c["gamma"] == -10)
def orbits(sites):
    d = {}
    for size in range(min(sites, 2)+1):
        for s in combinations(range(sites), size):
            for p in product((1,2), repeat=size):
                if sum(p)>2: continue
                v=[0]*sites
                for i,e in zip(s,p): v[i]=e
                v=tuple(v); d.setdefault(min(v,v[::-1]),set()).add(v)
    return [sorted(d[k]) for k in sorted(d)]
require(c["orbits"] == [[list(v) for v in o] for o in orbits(8)])
patterns=[q for q in product((-1,0,1),repeat=8) if sum(q)==0 and any(q)]
def feat(q,o):
    D=sum(x*x for x in q)//2
    return D*sum(prod(x**p for x,p in zip(q,pw)) for pw in o)
W=[]; K=[]; mult=[]
for q in patterns:
    D=sum(x*x for x in q)//2; v=[feat(q,o) for o in c["orbits"]]; penalty=[0]*len(v)
    for i in range(7):
        aa,bb=q[i]+1,q[i+1]+1
        targets=[]
        if abs(aa-bb)==1:
            t=list(q);t[i],t[i+1]=t[i+1],t[i];targets.append((tuple(t),1))
        elif {aa,bb}=={0,2}:
            t=list(q);t[i]=t[i+1]=0;targets.append((tuple(t),2))
        elif aa==bb==1:
            for z in (-1,1):
                t=list(q);t[i]=z;t[i+1]=-z;targets.append((tuple(t),1))
        for t,scale in targets:
            for j,o in enumerate(c["orbits"]): penalty[j]+=scale*feat(t,o)
    W.append(v); K.append([4*D*x-p for x,p in zip(v,penalty)]); mult.append(comb(8-2*D,(8-2*D)//2))
require(W == c["W"] and K == c["K"])
require([sum(mult[i]*W[i][j] for i in range(len(patterns))) for j in range(len(c["orbits"]))] == c["mean_num"])
require(c["mean_den"] == sum(mult))
h = json.loads((Path(__file__).parent.parent / "hamiltonian.json").read_text())
digest = hashlib.sha256(json.dumps(h, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
require(digest == "78365286838362d1b591e7afcef0f1e2a085fafc445e961f3b7e19603dbf1c6c", "Hamiltonian binding mismatch")
require(all(type(row) is int and 0 <= row < len(patterns) for row in c["support"]))
n = len(c["orbits"]); m = len(c["W"]); den = F(c["mean_den"])
y = [F(v) for v in c["y"]]; a = F(c["a"])
require(len(y) == len(c["support"]) and all(v >= 0 for v in y))
require(sum(y) == 1 and a < 0)
for j in range(n):
    lhs = sum(y[i] * (c["K"][row][j] + 10*c["W"][row][j]) for i,row in enumerate(c["support"]))
    require(lhs == a * F(c["mean_num"][j], den))
print(json.dumps({"verified": True, "support":len(y), "a":str(a), "scope":c["scope"]}))
