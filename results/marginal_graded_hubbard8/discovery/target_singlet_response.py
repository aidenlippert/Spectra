"""Fit the common response polynomial to the limiting Schur direction."""
from fractions import Fraction as F
from pathlib import Path
from math import lcm
import json,sys,time
import numpy as np
from scipy.linalg import eigh
from scipy.optimize import minimize
from numpy.polynomial import Polynomial,polynomial,chebyshev
sys.path.insert(0,str(Path(__file__).resolve().parents[3]))
from experiments.marginal_symmetry_moments import projected_moments
from results.marginal_graded_hubbard8.discovery.discover_singlet_moment_energy import chebyshev as integer_chebyshev
from results.marginal_graded_hubbard8.discovery.singlet_moment_energy import replay
ROOT=Path(__file__).resolve().parents[1]
def main():
    out=ROOT/'singlet_moment_targeted';out.mkdir(exist_ok=True);started=time.monotonic()
    c=json.loads((ROOT/'singlet_moment_energy/certificate.json').read_text());M=json.loads((ROOT/'symmetry_moments/moments.json').read_text());K=projected_moments(M)
    tau=F('-4.236');gamma=F('-3.81');delta=float(gamma-tau);P=integer_chebyshev(7)
    factors=[tau*gamma,-(tau+gamma),F(1)];den=lcm(*(x.denominator for x in factors));factors=[int(x*den) for x in factors]
    Dm=np.empty((8,8,14,14));Lm=np.empty((8,14,14))
    for i,p in enumerate(P):
        for r in range(14):
            for s in range(14):
                Lm[i,r,s]=sum(a*(100*K[k+1][r][s]+381*K[k][r][s]) for k,a in enumerate(p))/(100*12**i)
        for j,q in enumerate(P):
            for r in range(14):
                for s in range(14):
                    Dm[i,j,r,s]=sum(a*b*sum(factors[z]*K[k+l+z][r][s] for z in range(3)) for k,a in enumerate(p) for l,b in enumerate(q))/(den*12**(i+j))
    G=np.array(M[0],float);base=-float(tau)*G-np.array(K[0],float)/delta
    mono=np.array([float(F(x))/1e18 for x in c['response_polynomial']]);composed=Polynomial(mono)(Polynomial([8,12]));initial=chebyshev.poly2cheb(composed.coef);initial=initial/initial[0]
    history=[]
    def objective(x):
        a=np.r_[1.,x];D=np.einsum('i,j,ijab->ab',a,a,Dm);L=np.einsum('i,iab->ab',a,Lm)
        X=np.linalg.solve(D,L);S=base+L.T@X/delta;values,vectors=eigh(S,G);v=vectors[:,0];w=X@v
        gradient=np.array([(2*w@(Lm[i]@v)-2*w@(np.einsum('j,jab->ab',a,Dm[i])@w))/delta for i in range(8)])
        return -values[0],-gradient[1:]
    initial_value=objective(initial[1:])[0]
    def callback(x):history.append(float(-objective(x)[0]))
    result=minimize(objective,initial[1:],jac=True,method='L-BFGS-B',bounds=[(-10,10)]*7,callback=callback,options={'maxiter':200,'ftol':1e-15,'gtol':1e-10,'maxls':40})
    a=np.r_[1.,result.x]
    coefficients=[str(round(sum(F(float(w))*p[k]/12**i for i,(w,p) in enumerate(zip(a,P)) if k<len(p))*10**18)) for k in range(8)]
    diagnostic={'target':str(tau),'initial_minimum_schur':float(-initial_value),'final_minimum_schur':float(-result.fun),'optimizer_success':bool(result.success),'message':str(result.message),'iterations':int(result.nit),'history':history,'response_polynomial':coefficients}
    (out/'numerical.json').write_text(json.dumps(diagnostic,indent=2)+'\n');print(diagnostic,flush=True)
    c['response_polynomial']=coefficients;c['lower']=str(tau)
    (out/'candidate.json').write_text(json.dumps(c,indent=2)+'\n')
    receipt=replay(c);(out/'certificate.json').write_text(json.dumps(c,indent=2)+'\n');(out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print('EXACT ACCEPTED',receipt['lower'],receipt['upper_float'],'seconds',time.monotonic()-started,flush=True)
if __name__=='__main__':main()
