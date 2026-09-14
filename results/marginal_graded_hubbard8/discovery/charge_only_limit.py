#!/usr/bin/env python3
"""Enumerate the finite charge-only transition matrix and emit/replay PF bounds.

The verifier path is deliberately stdlib-only: it reconstructs every pattern and
transition into an integer sparse dictionary and checks the supplied witnesses.
"""
from itertools import product
from pathlib import Path
import hashlib, json, sys

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "results/marginal_graded_hubbard8/charge_only_limit"
HAMILTONIAN = ROOT / "results/marginal_graded_hubbard8/hamiltonian.json"
HAMILTONIAN_DIGEST = "78365286838362d1b591e7afcef0f1e2a085fafc445e961f3b7e19603dbf1c6c"

def patterns(sites=8):
    return [q for q in product((-1, 0, 1), repeat=sites) if sum(q) == 0 and any(q)]

def build(sites=8):
    ps = patterns(sites); ix = {q:i for i,q in enumerate(ps)}; A = {}
    for r,q in enumerate(ps):
        A[(r,r)] = -4 * (sum(x*x for x in q)//2)
        for i in range(sites-1):
            a,b = q[i]+1,q[i+1]+1; targets=[]
            if abs(a-b)==1:
                z=list(q); z[i],z[i+1]=z[i+1],z[i]; targets.append((tuple(z),1))
            elif {a,b}=={0,2}:
                z=list(q); z[i]=z[i+1]=0; targets.append((tuple(z),2))
            elif a==b==1:
                for direction in (-1,1):
                    z=list(q); z[i],z[i+1]=direction,-direction; targets.append((tuple(z),1))
            for z,w in targets:
                if z in ix: A[(r,ix[z])] = A.get((r,ix[z]),0) + w
    return ps,A

def export():
    import numpy as np
    from scipy.sparse.linalg import eigs
    ps,A=build(); n=len(ps); M=np.zeros((n,n))
    for (r,c),v in A.items(): M[r,c]=v
    lam=float(eigs(M,k=1,which='LR',return_eigenvectors=False)[0].real)
    _,vr=eigs(M,k=1,which='LR'); _,vl=eigs(M.T,k=1,which='LR')
    scale=10**16; w=np.maximum(1,np.rint(np.abs(vr[:,0].real)*scale).astype(object)); y=np.maximum(1,np.rint(np.abs(vl[:,0].real)*scale).astype(object))
    # Exact integer checks use numerators after multiplying gamma by 100.
    def right(g,wv): return [sum(-100*A.get((r,c),0)*wv[c] + (((287 if g==-287 else 286)*wv[r]) if r==c else 0) for c in range(n)) for r in range(n)]
    def left(g,yv): return [sum(yv[r]*(-100*A.get((r,c),0) + ((287 if g==-287 else 286) if r==c else 0)) for r in range(n)) for c in range(n)]
    digest=hashlib.sha256(SOURCE.read_bytes()).hexdigest()
    cert={'sites':8,'patterns':n,'nnz':len(A),'matrix_convention':'A=-4D+positive charge transition rows','hamiltonian_binding':{'model':'graded Hubbard','U':4,'t':1,'artifact':'results/marginal_graded_hubbard8/hamiltonian.json','canonical_sha256':HAMILTONIAN_DIGEST},'pf_eigenvalue':lam,'gamma_units':'hundredths','gamma_right_numerator':287,'gamma_left_numerator':286,'scale':scale,'right_weight':list(map(int,w)),'left_weight':list(map(int,y)),'right_residual_gamma_-2.87':right(-287,w),'left_residual_gamma_-2.86':left(-286,y),'source_generator_sha256':digest}
    cert['checks']={'right_strict_positive':min(cert['right_residual_gamma_-2.87'])>0,'left_strict_negative':max(cert['left_residual_gamma_-2.86'])<0}
    OUT.mkdir(parents=True,exist_ok=True); (OUT/'certificate.json').write_text(json.dumps(cert,indent=2)+'\n')
    (OUT/'receipt.json').write_text(json.dumps({'status':'exported','certificate':'certificate.json','hamiltonian_canonical_sha256':HAMILTONIAN_DIGEST},indent=2)+'\n')

def replay():
    cert=json.loads((OUT/'certificate.json').read_text())
    if cert.get('sites') != 8 or cert.get('patterns') != 1106 or cert.get('gamma_units') != 'hundredths' or cert.get('gamma_right_numerator') != 287 or cert.get('gamma_left_numerator') != 286: raise ValueError('certificate metadata mismatch')
    binding=cert.get('hamiltonian_binding',{})
    canonical=json.dumps(json.loads(HAMILTONIAN.read_text()),sort_keys=True,separators=(',',':'))
    if binding.get('U') != 4 or binding.get('t') != 1 or binding.get('canonical_sha256') != HAMILTONIAN_DIGEST or hashlib.sha256(canonical.encode()).hexdigest() != HAMILTONIAN_DIGEST: raise ValueError('Hamiltonian binding mismatch')
    ps,A=build(8); n=len(ps)
    w,y=cert.get('right_weight'),cert.get('left_weight')
    for name,v in (('right_weight',w),('left_weight',y)):
        if not isinstance(v,list) or len(v)!=n or any(type(x) is not int or x<=0 for x in v): raise ValueError(name+' must contain exactly 1106 positive integers')
    def right(w): return [sum(-100*A.get((r,c),0)*w[c] + (287*w[r] if r==c else 0) for c in range(n)) for r in range(n)]
    def left(y): return [sum(y[r]*(-100*A.get((r,c),0) + (286 if r==c else 0)) for r in range(n)) for c in range(n)]
    rp,ln=right(w),left(y)
    if min(rp)<=0 or max(ln)>=0: raise ValueError('PF residual sign gate failed')
    result={'status':'verified','stdlib_only':True,'patterns':n,'nnz':len(A),'right_min':min(rp),'left_max':max(ln),'right_strict_positive':True,'left_strict_negative':True,'hamiltonian_canonical_sha256':HAMILTONIAN_DIGEST}
    (OUT/'independent_replay.json').write_text(json.dumps(result,indent=2)+'\n'); print(json.dumps(result,indent=2))

if __name__=='__main__':
    if '--replay' in sys.argv: replay()
    else: export()
