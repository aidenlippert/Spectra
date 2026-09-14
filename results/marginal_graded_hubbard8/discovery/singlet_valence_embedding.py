"""Exact 14-dimensional noncrossing valence-bond embedding for H8.

This is a bounded geometry artifact.  It does not alter production energy gates.
All arithmetic used for the certificate and replay is integer/rational stdlib.
"""
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import argparse, hashlib, json
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from experiments.marginal_determinant_tree import DeterminantOracle
from experiments.marginal_transfer_verify import apply_word

ROOT = Path(__file__).resolve().parents[1]

def matchings(lo=0, hi=8):
    if lo == hi: return [()]
    out=[]
    for j in range(lo+1, hi, 2):
        for a in matchings(lo+1,j):
            for b in matchings(j+1,hi): out.append(((lo,j),)+a+b)
    return out

def valence_states():
    return [sum(1<<(2*i) for i in up)|sum(1<<(2*i+1) for i in range(8) if i not in up)
            for up in combinations(range(8),4)]

def singlet_vector(pairs):
    # Each pair contributes up_i down_j - down_i up_j; creation operators are
    # canonicalized by sorting occupied modes, introducing the CAR sign.
    v={}
    for choices in __import__('itertools').product((0,1), repeat=4):
        modes=[]; sign=1
        for (i,j), c in zip(pairs,choices):
            modes += [2*i+(c==1), 2*j+(c==0)]
            if c: sign=-sign
        inv=sum(modes[a]>modes[b] for a in range(8) for b in range(a+1,8))
        s=sign*(-1 if inv%2 else 1); state=sum(1<<m for m in modes)
        v[state]=v.get(state,0)+s
    return {s:c for s,c in v.items() if c}

def gram(V, basis): return [[sum(V[k].get(i,0)*V[l].get(i,0) for i in basis) for l in range(len(V))] for k in range(len(V))]
def rank(a):
    a=[[F(x) for x in r] for r in a]; m=len(a); n=len(a[0]) if m else 0; r=0
    for c in range(n):
        p=next((i for i in range(r,m) if a[i][c]),None)
        if p is None: continue
        a[r],a[p]=a[p],a[r]; q=a[r][c]; a[r]=[x/q for x in a[r]]
        for i in range(m):
            if i!=r and a[i][c]:
                q=a[i][c]; a[i]=[x-q*y for x,y in zip(a[i],a[r])]
        r+=1
    return r
def ldl(a):
    n=len(a); d=[]; L=[[F(int(i==j)) for j in range(n)] for i in range(n)]
    for i in range(n):
        x=F(a[i][i])-sum(L[i][k]*L[i][k]*d[k] for k in range(i)); d.append(x)
        if not x: raise ValueError('singular Gram/Schur matrix')
        for j in range(i+1,n): L[j][i]=(F(a[j][i])-sum(L[j][k]*L[i][k]*d[k] for k in range(i)))/x
    return d
def jsonable(x):
    if isinstance(x,F): return str(x)
    if isinstance(x,dict): return {str(k):jsonable(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)): return [jsonable(v) for v in x]
    return x

def main(out):
    out=Path(out); out.mkdir(parents=True,exist_ok=True)
    h=json.loads((ROOT/'hamiltonian.json').read_text()); h['modes']=16; h['particles']=8
    digest=hashlib.sha256(json.dumps(h,sort_keys=True,separators=(',',':')).encode()).hexdigest()
    p=valence_states(); ms=matchings(); V=[singlet_vector(x) for x in ms]
    # S+ maps each down mode to up mode, with exact CAR signs.
    sz1=[sum(1<<(2*i) for i in up)|sum(1<<(2*i+1) for i in dn) for up in combinations(range(8),5) for dn in combinations(range(8),3) if not set(up)&set(dn)]
    splus=[]
    for s in p:
        row={}
        for i in range(8):
            z=apply_word(((1,2*i),(0,2*i+1)),s)
            if z: row[z[0]]=row.get(z[0],0)+z[1]
        splus.append(row)
    sr=[[row.get(t,0) for t in sz1] for row in splus]
    G=gram(V,p)
    oracle=DeterminantOracle(h); cols=[{t:x for t,x in oracle.action(s).items() if t not in set(p)} for s in p]
    W=[[sum(a.get(t,0)*b.get(t,0) for t in cols[i]) for b in cols] for i,a in enumerate(cols)]
    WG=[[sum(V[i].get(p[k],0)*W[k][l]*V[j].get(p[l],0) for k in range(70) for l in range(70)) for j in range(14)] for i in range(14)]
    tau=F('-21/5'); gamma=F('-381/100'); den=gamma-tau
    S=[[-tau*G[i][j]-WG[i][j]/den for j in range(14)] for i in range(14)]
    payload={'schema':'h8-singlet-valence-embedding-v1','sites':8,'matching_count':len(ms),'matchings':ms,'basis':V,'valence_states':p,'gram':G,'gram_ldl':[str(x) for x in ldl(G)],'splus_rank':rank(sr),'splus_kernel_dimension':70-rank(sr),'splus_target_dimension':len(sz1),'leakage_gram':W,'projected_leakage':WG,'gamma':str(gamma),'tau':str(tau),'restricted_schur_ldl':[str(x) for x in ldl(S)],'restricted_schur_negative_pivots':sum(x<0 for x in ldl(S)),'hamiltonian_sha256':digest,'scope':'Finite 14-dimensional noncrossing valence-bond geometry; geometry diagnostic only, no energy interval certificate.','model_assumptions':'U=4, t=1 real nearest-neighbor bipartite chain A4/B4; positive connected bipartite half-filled Hubbard qualifies for Lieb theorem (PhysRevLett.62.1201), recorded as an assumption only.'}
    (out/'certificate.json').write_text(json.dumps(jsonable(payload),indent=2)+'\n')
    # The independently authored read-only verifier is maintained separately.
    (out/'source_h_digest.json').write_text(json.dumps({'sha256':digest},indent=2)+'\n')
if __name__=='__main__':
    ap=argparse.ArgumentParser(); ap.add_argument('--out',default=str(ROOT/'singlet_valence_embedding')); main(ap.parse_args().out)
