"""Numerical proposals for three explicit extended-Hubbard targets."""
from pathlib import Path
from fractions import Fraction as F
import json,sys,math,time
import numpy as np
from scipy.linalg import eigh
from scipy.optimize import minimize
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT))
from experiments.marginal_local_hubbard_block import sector_matrices
from experiments.marginal_boundary_transfer import compile_block
BASE=ROOT/'results/marginal_graded_hubbard8'

def main():
    source=json.loads((BASE/'tiled_eight_upper/certificate.json').read_text())
    window=source['lower_window'];results=[];start=time.monotonic()
    for U,t,V in [(F(4),F(1),F(1,2)),(F(4),F(1),F(-1,2)),(F(3),F(2,3),F(1,2))]:
        onsite=[F(x)*U/4 for x in window['onsite_profile']]
        hopping=[F(x)*t for x in window['hopping_profile']]
        density=[F(x)*V for x in window['hopping_profile']]
        sectors=sector_matrices(6,5*U/6,t,onsite,hopping,V,density)
        lows=[]
        for key,(_,matrix,columns) in sectors.items():
            weights=np.sqrt([sum(a*a for a in v.values()) for v in columns])
            A=np.array(matrix,dtype=float)/weights[:,None]/weights[None,:]
            lows.append((float(eigh(A,subset_by_index=[0,0],eigvals_only=True)[0]),list(key)))
        low,key=min(lows)
        lower=F(math.floor(low*100000)-1,100000)
        certificate={'kind':'local_hubbard_block_v1','sites':6,'U':str(5*U/6),'t':str(t),'V':str(V),
            'onsite_profile':list(map(str,onsite)),'hopping_profile':list(map(str,hopping)),
            'density_profile':list(map(str,density)),'lower':str(lower)}
        c=compile_block(source['upper'],source['hamiltonian'],U,t,V);norm=c['block_norm']
        G=np.array([[float(F(x,norm)) for x in row] for row in c['G']])
        J=np.array([[float(F(x,norm)) for x in row] for row in c['J']])
        kernels=[(m,n,np.array(B,dtype=float),np.array([[float(F(x,norm*norm)) for x in row] for row in W])) for m,n,B,W in c['kernels']]
        def energy(x):
            z=[1,*x];B=sum(z[m]*z[n]*B for m,n,B,W in kernels);W=sum(z[m]*z[n]*W for m,n,B,W in kernels)
            T=B@G;values,right=np.linalg.eig(T);idx=np.argmax(values.real);lam=values[idx]
            vals,left=np.linalg.eig(T.T);l=left[:,np.argmin(abs(vals-lam))];r=right[:,idx]
            answer=((l@(B@J)@r)/lam+(l@(B@W)@r)/(lam*lam))/(l@r)/8
            if abs(answer.imag)>1e-9 or lam.real<=0:raise ValueError('Unsupported numerical transfer spectrum')
            return float(answer.real)
        opt=minimize(energy,[.356039,.177774],method='Nelder-Mead',options={'maxiter':250,'xatol':1e-10,'fatol':1e-12})
        a,b=[format(float(x),'.6f') for x in opt.x]
        result={'target':{'U':str(U),'t':str(t),'V':str(V)},'a':a,'b':b,'lower_certificate':certificate,
            'numerical_local_minimum':low,'numerical_active_sector':key,
            'numerical_upper_density':energy([float(a),float(b)]),'optimizer_success':bool(opt.success)}
        results.append(result);print(json.dumps(result),flush=True)
    out=BASE/'density_transfer';out.mkdir(exist_ok=True)
    (out/'numeric_proposals.json').write_text(json.dumps({'targets':results,'seconds':time.monotonic()-start,
        'scope':'Numerical proposal only. Lower profiles reused and scaled, not optimized for each target; all-Fock PSD replay required.'},indent=2)+'\n')

if __name__=='__main__':main()
