"""Working-copy adapters; supplied project files are never modified."""
from pathlib import Path
import sys, json, time, os
ROOT=Path(os.environ.get('SPECTRA_ROOT', '/mnt/data/spectra_work/spectra_h8')).resolve()
OUT=Path(os.environ.get('SPECTRA_OUTPUT', str(Path(__file__).resolve().parent/'runs'))).resolve()
OUT.mkdir(parents=True,exist_ok=True)
sys.path.insert(0,str(ROOT))
import numpy as np
from scipy import linalg,sparse
from scipy.sparse.linalg import LinearOperator,cg
from research.sector_quotient_20260914.search import Operator,export
from research.sector_quotient_20260914.eliminated import Quotient
from research.sector_quotient_20260914.boundary import positive
PREP=ROOT/'results/sector_quotient_20260914/prepared_linear_closure'
CACHE=ROOT/'results/sector_quotient_20260914/normal_linear_closure'
SOURCE=ROOT/'results/sector_quotient_20260914/candidates/optimize_number_linear_closure'
U=-9.254878154398348813867703762444199296

def load_operator():
    meta=json.loads((PREP/'frame.json').read_text())
    def rebase(v):
        if isinstance(v,str) and v.startswith('/Users/aidenlippert/Documents/Spectra/'):
            return str(ROOT/v.split('/Users/aidenlippert/Documents/Spectra/',1)[1])
        if isinstance(v,list): return [rebase(x) for x in v]
        if isinstance(v,dict): return {k:rebase(x) for k,x in v.items()}
        return v
    return Operator(PREP,rebase(meta))

def load_source(op,source=SOURCE/'checkpoint.npz'):
    raw=np.load(source)
    Q=[raw[f'Q_{i}'].copy() for i in range(len(op.Q))]
    return Q,raw['y'].copy(),raw['x'].copy()

def make_inverse(quo,cache=CACHE):
    C=np.load(cache/'cholesky.npy')
    inverse=lambda v:linalg.cho_solve((C,True),v,check_finite=False)
    gv=inverse(quo.v);den=float(quo.v@gv)
    def conditioned(v):
        v=quo.project(v)
        a=inverse(v)-gv*(gv@v)/den
        return quo.project(a)
    return conditioned,C

def record_json(p,d):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(d,indent=2))

class FastResidual:
    """Floating proposal-only reconstruction of the full spin-averaged residual."""
    def __init__(self,op):
        from itertools import combinations
        from math import comb
        from scipy.sparse.linalg import splu
        self.op=op;self.comb=comb
        T=sparse.load_npz(PREP/'twirl.npz');sel=np.load(PREP/'selected.npy')
        self.Tc=T[:,sel].tocsc();self.factor=splu(T[sel][:,sel].tocsc());self.sel=sel;self.rows=op.meta['rows'];self.parts=[]
        for k in (1,2,3):
            inds=list(combinations(range(16),k));lookup={v:i for i,v in enumerate(inds)};par=list(range(len(inds)))
            def root(i):
                while par[i]!=i:par[i]=par[par[i]];i=par[i]
                return i
            triples=[]
            for idx,w in enumerate(self.rows):
                if len(w)!=2*k:continue
                i=lookup[tuple(p for c,p in w[:k])];j=lookup[tuple(p for c,p in w[k:])];triples.append((idx,i,j))
                a,b=root(i),root(j)
                if a!=b:par[a]=b
            groups={}
            for i in range(len(inds)):groups.setdefault(root(i),[]).append(i)
            components=[]
            for gg in groups.values():
                loc={i:j for j,i in enumerate(gg)};pairs=[(idx,loc[i],loc[j]) for idx,i,j in triples if i in loc]
                components.append((len(gg),np.array([q[0] for q in pairs],dtype=int),np.array([q[1] for q in pairs],dtype=int),np.array([q[2] for q in pairs],dtype=int)))
            self.parts.append((k,components))
    def evaluate(self,selected_residual):
        full=-self.Tc@self.factor.solve(np.asarray(selected_residual))
        if np.max(abs(full[self.sel]+selected_residual))>1e-7:raise AssertionError('Residual interpolation mismatch')
        result=float(full[0]);bounds=[]
        for k,components in self.parts:
            low=float('inf');l1=0.
            for n,ind,i,j in components:
                vals=((-1)**(k*(k-1)//2))*full[ind];M=np.zeros((n,n));M[i,j]=vals;M[j,i]=vals
                ev=float(linalg.eigvalsh(M,subset_by_index=(0,0),check_finite=False)[0]);low=min(low,ev);l1+=float(np.sum(abs(vals)*(1+(i!=j))))
            chosen=max(-l1,self.comb(8,k)*low);result+=chosen;bounds.append(chosen)
        return result,bounds
