"""Numerical singlet Krylov response probe (no energy certificate claim)."""
import json, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from experiments.marginal_determinant_tree import DeterminantOracle

ROOT=Path(__file__).resolve().parents[1]

def add(a,b,scale=1):
    c=dict(a)
    for s,x in b.items(): c[s]=c.get(s,0)+scale*x
    return {s:x for s,x in c.items() if x}
def act(o,v):
    z={}
    for s,a in v.items():
        for t,x in o.action(s).items(): z[t]=z.get(t,0)+a*x
    return {s:x for s,x in z.items() if x}
def qact(o,v,retained):
    return {s:x for s,x in act(o,v).items() if s not in retained}
def dot(a,b): return sum(x*b.get(s,0) for s,x in a.items())
def main():
    src=ROOT/'singlet_valence_embedding'/'certificate.json'; c=json.loads(src.read_text())
    h=json.loads((ROOT/'hamiltonian.json').read_text()); h.update(modes=16,particles=8)
    o=DeterminantOracle(h); p=c['valence_states']; V=[]
    for v in c['basis']: V.append({int(s):int(a) for s,a in v.items()})
    # H maps this valence space entirely into Q.
    retained=set(p); W=[qact(o,v,retained) for v in V]
    G=np.array([[dot(a,b) for b in V] for a in V],float)
    out=[]
    for order in range(0,6):
        B=[]
        for v in V:
            x=dict(W[V.index(v)]); y={}
            alpha=20.; den=alpha-(-4.3)
            for k in range(order+1):
                y=add(y,x,1./den)
                x=add({s:alpha*a for s,a in x.items()},qact(o,x,retained),-1.)
                x={s:a/den for s,a in x.items()}
            B.append(y)
        # remove numerical linear dependence by Gram eigenspectrum
        GB=np.array([[dot(a,b) for b in B] for a in B],float)
        ev=np.linalg.eigvalsh(GB); keep=np.where(ev>1e-8*max(ev))[0]
        # use eigenvectors as an orthonormal coordinate system in B
        vals,vecs=np.linalg.eigh(GB); T=vecs[:,vals>1e-8*max(vals)]/np.sqrt(vals[vals>1e-8*max(vals)])
        def comb(col):
            z={}
            for a,w in zip(B,col): z=add(z,a,float(w))
            return z
        E=[comb(T[:,j]) for j in range(T.shape[1])]
        try:
            HE=[qact(o,x,retained) for x in E]
        except ValueError as exc:
            out.append({'order':order,'raw_directions':len(B),'independent_directions':len(E),'blocker':str(exc),'cached_actions':len(o.cache)})
            break
        metric=np.array([[dot(a,b) for b in E] for a in E],float)
        hbb=np.array([[dot(a,b) for b in HE] for a in E],float)
        squared=np.array([[dot(a,b) for b in HE] for a in HE],float)
        C=np.array([[dot(a,b) for b in W] for a in E],float)
        HC=np.array([[dot(a,b) for b in W] for a in HE],float)
        WW=np.array([[dot(a,b) for b in W] for a in W],float)
        tau=-4.3; gamma=-3.81; delta=gamma-tau
        D=squared-(tau+gamma)*hbb+tau*gamma*metric
        # Formula equivalent to current generalized response with G=V^TV.
        L=HC-gamma*C
        corr=L.T@np.linalg.solve(D,L)
        S=-tau*G-(WW-corr)/delta
        Lg=np.linalg.cholesky(G); se=np.linalg.eigvalsh(np.linalg.solve(Lg,S) @ np.linalg.inv(Lg.T))
        out.append({'order':order,'raw_directions':len(B),'independent_directions':len(E),'min_schur_eigenvalue':float(se.min()),'negative_count':int((se< -1e-8).sum()),'max_support':max(map(len,E)),'response_support_total':sum(map(len,E))})
    dest=ROOT/'singlet_response_discovery'; dest.mkdir(exist_ok=True)
    (dest/'diagnostics.json').write_text(json.dumps({'schema':'h8-singlet-krylov-response-diagnostic-v1','source':str(src),'gamma':'-3.81','tau':'-4.3','alpha':'20','polynomial':'Neumann p_d(A)=sum_{k=0}^d ((alpha-A)/(alpha-tau))^k/(alpha-tau)','results':out,'scope':'Numerical response probe only; no exact certificate or energy interval.'},indent=2)+'\n')
if __name__=='__main__': main()
