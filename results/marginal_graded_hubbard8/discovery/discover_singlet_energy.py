"""Fit a common short inverse polynomial; acceptance is separate exact replay."""
from pathlib import Path
from fractions import Fraction as F
import json,sys,copy,math
import numpy as np
from numpy.polynomial import Polynomial as Poly
from scipy.linalg import eigh
sys.path.insert(0,str(Path(__file__).resolve().parents[3]))
from results.marginal_graded_hubbard8.discovery.singlet_energy_replay import prepare,apply,gram,replay
ROOT=Path(__file__).resolve().parents[1]
def combine(vectors,coefficients):
    out={}
    for v,a in zip(vectors,coefficients):
        for s,x in v.items():out[s]=out.get(s,0)+a*x
    return {s:x for s,x in out.items() if x}
def numeric_schur(G,W,B,AB,gamma,tau):
    arr=lambda a,b:np.array(gram(a,b),float)
    D=arr(AB,AB)-(tau+gamma)*arr(B,AB)+tau*gamma*arr(B,B)
    L=arr(AB,W)-gamma*arr(B,W)
    S=-tau*G-(arr(W,W)-L.T@np.linalg.solve(D,L))/(gamma-tau)
    return eigh(S,G,eigvals_only=True)[0]
def main():
    out=ROOT/'singlet_energy_fit';out.mkdir(exist_ok=True)
    c=json.loads((ROOT/'singlet_energy_first/certificate.json').read_text())
    data,o,V,W,_,_=prepare(c);P=set(c['embedding']['valence_states']);G=np.array(data['G'],float)
    # Stable Chebyshev coordinates, fit is a proposal and assumes no spectral proof.
    center=8;scale=12;polys=[Poly([1]),Poly([-center/scale,1/scale])]
    vectors=[W];images=[[apply(o,w,P) for w in W]];records=[];best=None
    for degree in range(6):
        if degree:
            if degree==1:
                layer=[combine([a,w],[F(1,scale),F(-center,scale)]) for a,w in zip(images[0],W)]
            else:
                layer=[combine([a,w,old],[F(2,scale),F(-2*center,scale),-1]) for a,w,old in zip(images[-1],vectors[-1],vectors[-2])]
                polys.append(2*polys[1]*polys[-1]-polys[-2])
            vectors.append(layer);images.append([apply(o,w,P) for w in layer])
        # Minimize sum of squared residuals over the 14 independent RHS columns.
        tau=-4.3
        R=[[combine([a,w],[1.,-tau]) for a,w in zip(aa,ww)] for aa,ww in zip(images,vectors)]
        M=np.array([[sum(sum(float(x)*float(b.get(s,0)) for s,x in a.items()) for a,b in zip(left,right)) for right in R] for left in R])
        rhs=np.array([sum(sum(float(x)*float(w.get(s,0)) for s,x in a.items()) for a,w in zip(left,W)) for left in R])
        coeff=np.linalg.solve(M,rhs)
        B=[combine([v[j] for v in vectors],coeff) for j in range(14)]
        AB=[combine([v[j] for v in images],coeff) for j in range(14)]
        thresholds=[-4.3,-4.4,-4.5,-4.6,-4.8,-5.,-5.2,-5.5,-6.]
        scan=[{'tau':t,'min_schur':float(numeric_schur(G,W,B,AB,-3.81,t))} for t in thresholds]
        poly=sum((p*a for p,a in zip(polys,coeff)),Poly([0]))
        rational=[str(F(float(x)).limit_denominator(10**10)) for x in poly.coef]
        row={'degree':degree,'chebyshev_coefficients':coeff.tolist(),'polynomial':rational,'scan':scan,'unique_action_states':len(o.cache),'referenced_determinants':o.referenced_state_count()}
        records.append(row);print(row,flush=True)
        accepted=next((s['tau'] for s in scan if s['min_schur']>1e-5),None)
        if accepted is not None and (best is None or accepted>best[0]):best=(accepted,rational,degree)
        (out/'numerical.json').write_text(json.dumps(records,indent=2)+'\n')
    require_best=best is not None
    if not require_best:raise ValueError('No numerically feasible polynomial')
    c['lower']=str(best[0]);c['response_polynomial']=best[1]
    c['independent_upper']=json.loads((ROOT/'singlet_upper_krylov/upper_degree4.json').read_text())
    (out/'candidate.json').write_text(json.dumps(c,indent=2)+'\n')
    receipt=replay(c)
    (out/'certificate.json').write_text(json.dumps(c,indent=2)+'\n');(out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print('EXACT ACCEPTED',receipt['lower'],receipt['upper_float'],flush=True)
if __name__=='__main__':main()
