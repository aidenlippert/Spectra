"""Small charge-conserving MPS variational search on frozen CAR fixtures.

The tensor parameterization enforces particle number at every amplitude:
virtual charge q advances by physical occupation n.  The final exact witness
export is intentionally exponential and is reported separately.
"""
import sys,json,hashlib,time
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
ROOT=Path(__file__).resolve().parents[3]; sys.path.insert(0,str(ROOT))
from experiments.marginal_determinant_tree import DeterminantOracle

def charge_shapes(L,N,D):
    # each cut has charge sectors 0..N, capped by D; deterministic dimensions
    return [min(D,q+1,N-q+1) for q in range(N+1)]
def unpack(x,L,N,D):
    ds=charge_shapes(L,N,D); out=[]; k=0
    for i in range(L):
        qmax=min(i,N); rmax=min(i+1,N); blocks=[]
        for q in range(qmax+1):
            for n in (0,1):
                rr=q+n
                if rr>rmax: continue
                a=ds[q]; b=ds[rr]; z=a*b
                blocks.append((q,n, x[k:k+z].reshape(a,b))); k+=z
        out.append(blocks)
    return out,k
def amplitudes(x,L,N,D):
    ts,_=unpack(x,L,N,D); vec=np.zeros(1<<L)
    for s in range(1<<L):
        if s.bit_count()!=N: continue
        v=np.array([1.]) ; q=0
        for i in range(L):
            n=(s>>i)&1; found=[a for a in ts[i] if a[0]==q and a[1]==n]
            if not found: v=None; break
            v=v@found[0][2]; q+=n
        if v is not None: vec[s]=v[0]
    return vec
def term_expect(x,word,L,N,D):
    """Direct MPS transfer contraction of a CAR word using JW local factors."""
    ts,_=unpack(x,L,N,D); I=np.eye(2); Z=np.diag([1.,-1.]); cp=np.array([[0.,0.],[1.,0.]]); an=cp.T
    fac=[I.copy() for _ in range(L)]
    # word is operator product left-to-right; JW parity is on lower modes.
    for c,j in word:
        op=cp if c else an
        for k in range(j): fac[k]=fac[k]@Z
        fac[j]=fac[j]@op
    E=np.array([[1.]])
    for i in range(L):
        ds=charge_shapes(L,N,D); ql=list(range(min(i,N)+1)); qr=list(range(min(i+1,N)+1)); dl=sum(ds[q] for q in ql); dr=sum(ds[q] for q in qr)
        # recover padded physical blocks from charge-sector tensors
        A=[np.zeros((dl,dr)) for _ in (0,1)]; lo={q:sum(ds[z] for z in ql if z<q) for q in ql}; ro={q:sum(ds[z] for z in qr if z<q) for q in qr}
        for q,n,b in ts[i]: A[n][lo[q]:lo[q]+b.shape[0],ro[q+n]:ro[q+n]+b.shape[1]]=b
        E=sum(A[s].T@E@A[t]*fac[i][s,t] for s in (0,1) for t in (0,1))
    ds=charge_shapes(L,N,D); off=sum(ds[q] for q in range(N)); return E[off,off]
def contracted_energy(x,terms,L,N,D):
    norm=term_expect(x,(),L,N,D)
    return sum(float(c)*term_expect(x,w,L,N,D) for w,c in terms)/norm
def load(name):
    p=ROOT/f'results/certificate_scaling/active_space_ladder/{name}/fixture.json'; return p,json.loads(p.read_text())
def run(name,D=2):
    path,fix=load(name); oracle=DeterminantOracle(fix); L=fix['modes']; N=fix['particles']; states=[s for s in range(1<<L) if s.bit_count()==N]; ix={s:i for i,s in enumerate(states)}
    H=np.zeros((len(states),len(states))); terms=list(oracle.h.items())
    for s in states:
        for t,c in oracle.action(s).items(): H[ix[t],ix[s]]=float(c)
    rng=np.random.default_rng(20260913); _,n=unpack(np.empty(100000),L,N,D); x=rng.normal(0,.1,n)
    def f(y): return contracted_energy(y,terms,L,N,D)
    t=time.time(); r=minimize(f,x,method='BFGS',options={'maxiter':80,'gtol':1e-7}); w=amplitudes(r.x,L,N,D)[states]; w/=np.linalg.norm(w)
    den=10**10; a=np.rint(w*den).astype(np.int64); nz=np.flatnonzero(a); wit={'states':[states[i] for i in nz],'amplitudes':[int(a[s]) for s in nz]}; ex=oracle.upper(wit)
    return {'fixture':str(path),'fixture_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'modes':L,'particles':N,'bond_dimension':D,'parameter_count':n,'wall_seconds':time.time()-t,'optimizer_success':bool(r.success),'message':r.message,'float_energy':float(f(r.x)),'exact_integer_witness_energy':str(ex),'support':len(nz),'export_is_exponential':True,'method':'charge-conserving MPS; H matrix used only for small-fixture optimization; witness replay exact CAR'} ,wit
def main():
    out=ROOT/'results/wave2_20260913/tensor_network'; out.mkdir(parents=True,exist_ok=True); rows=[]
    for name in ('h4','h6'):
        try:
            row,w=run(name); rows.append(row); (out/f'{name}_witness.json').write_text(json.dumps(w,indent=2)+'\n')
        except Exception as e: rows.append({'fixture':name,'status':'failed','error':repr(e)})
    (out/'report.json').write_text(json.dumps({'rows':rows,'scope':'Wave2 charge-conserving MPS benchmark'},indent=2)+'\n'); print(json.dumps(rows,indent=2))
if __name__=='__main__': main()
