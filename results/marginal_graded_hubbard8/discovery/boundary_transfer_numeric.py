"""Finite H8 numerical gate proposal from exact local transfer tensors."""
from pathlib import Path
import sys,json,time
import numpy as np
from scipy.optimize import minimize
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT))
from experiments.marginal_boundary_transfer import compile_block,contract

def main():
    source=json.loads((ROOT/'results/marginal_graded_hubbard8/tiled_eight_upper/certificate.json').read_text())
    start=time.monotonic();c=compile_block(source['upper'],source['hamiltonian'])
    norm=c['block_norm']
    G=np.array(c['G'],dtype=float)/float(norm)
    J=np.array(c['J'],dtype=float)/float(norm)
    kernels=[(m,n,np.array(B,dtype=float),np.array(W,dtype=float)/float(norm*norm)) for m,n,B,W in c['kernels']]
    def energy(x):
        a,b=x;z=[1,a,b]
        B=sum(z[m]*z[n]*B for m,n,B,W in kernels)
        W=sum(z[m]*z[n]*W for m,n,B,W in kernels)
        T=B@G;values,right=np.linalg.eig(T);idx=np.argmax(values.real);lam=values[idx]
        lv,left=np.linalg.eig(T.T);li=np.argmin(abs(lv-lam))
        r=right[:,idx];l=left[:,li];overlap=l@r
        answer=((l@(B@J)@r)/lam+(l@(B@W)@r)/(lam*lam))/overlap/8
        if abs(answer.imag)>1e-9 or lam.real<=0:raise ValueError('Numerical transfer spectrum unsupported')
        return float(answer.real)
    opt=minimize(energy,[.229108,.02],method='Nelder-Mead',options={'maxiter':250,'xatol':1e-10,'fatol':1e-12})
    a,b=[format(float(x),'.6f') for x in opt.x]
    numeric={'a':a,'b':b,'numerical_density':energy([float(a),float(b)]),
             'optimizer_success':bool(opt.success),'scope':'Floating-point thermodynamic proposal only. Fixed16-dimensional transfer spectrum; no certified large-chain energy.'}
    out=ROOT/'results/marginal_graded_hubbard8/boundary_transfer';out.mkdir(exist_ok=True)
    (out/'numeric_proposal.json').write_text(json.dumps(numeric,indent=2)+'\n')
    print(json.dumps(numeric),flush=True)
    for q in (1,2,3,8):
        r=contract(c,a,b,q)
        print(q,float(__import__('fractions').Fraction(r['upper_per_site'])),flush=True)
    print('seconds',time.monotonic()-start,flush=True)

if __name__=='__main__':main()
