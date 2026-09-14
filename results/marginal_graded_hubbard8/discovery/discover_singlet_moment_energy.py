"""Propose compact response and variational polynomials using exact moments.

Saved moment tables are proposal inputs only. Exact acceptance independently
recomputes moments from the original Hamiltonian via singlet_moment_energy.
"""
from pathlib import Path
from fractions import Fraction as F
import json,sys,time
import numpy as np
from scipy.linalg import eigh
sys.path.insert(0,str(Path(__file__).resolve().parents[3]))
from results.marginal_graded_hubbard8.discovery.singlet_moment_energy import response_data,polynomial_upper,upper_monomials,replay
from experiments.marginal_symmetry_moments import projected_moments
ROOT=Path(__file__).resolve().parents[1]
def add(a,b):
    return [(a[k] if k<len(a) else 0)+(b[k] if k<len(b) else 0) for k in range(max(len(a),len(b)))]
def times_x(a):return [0]+a

def chebyshev(maximum):
    P=[[1],[-8,1]]
    for k in range(2,maximum+1):P.append(add(add([2*x for x in times_x(P[-1])],[-16*x for x in P[-1]]),[-144*x for x in P[-2]]))
    return P[:maximum+1]
def main(moment_source='symmetry_moments',output='singlet_moment_energy',response_degrees=(5,6,7),upper_degrees=(4,6,8),fitting_tau=F(-17,4),lower_grid=None,chebyshev_upper=False):
    started=time.monotonic();out=ROOT/output;out.mkdir(exist_ok=True)
    M=json.loads((ROOT/moment_source/'moments.json').read_text());K=projected_moments(M)
    trace=[sum(m[i][i] for i in range(14)) for m in K];P=chebyshev(max(*response_degrees,*upper_degrees))
    records=[];best=None
    for degree in response_degrees:
        # Exact integer contractions form the least-squares normal equations.
        tau=F(fitting_tau);d=tau.denominator;R=[add([d*x for x in times_x(p)],[-tau.numerator*x for x in p]) for p in P[:degree+1]]
        mat=np.array([[sum(a*b*trace[k+l] for k,a in enumerate(r) for l,b in enumerate(s))/(d*d*12**(i+j)) for j,s in enumerate(R)] for i,r in enumerate(R)])
        rhs=np.array([sum(a*trace[k] for k,a in enumerate(r))/(d*12**i) for i,r in enumerate(R)])
        weights=np.linalg.solve(mat,rhs)
        coeff=[sum(F(float(w))*p[k]/12**i for i,(w,p) in enumerate(zip(weights,P)) if k<len(p)) for k in range(degree+1)]
        integer=[str(round(x*10**18)) for x in coeff]
        data=response_data(M,integer,F(-381,100));arr=lambda key:np.array(data[key],float)
        G=arr('G');scan=[]
        for lower in (lower_grid or [-4.24,-4.25,-4.26,-4.28,-4.3,-4.4]):
            gamma=-3.81;D=arr('squared')-(lower+gamma)*arr('hamiltonian')+lower*gamma*arr('metric')
            L=arr('h_coupling')-gamma*arr('coupling');S=-lower*G-(arr('leakage')-L.T@np.linalg.solve(D,L))/(gamma-lower)
            minimum=float(eigh(S,G,eigvals_only=True)[0]);scan.append({'lower':lower,'minimum_schur':minimum})
        feasible=next((x['lower'] for x in scan if x['minimum_schur']>1e-5),None)
        row={'degree':degree,'integer_coefficients':integer,'scan':scan};records.append(row)
        if feasible is not None and (best is None or feasible>best[0]):best=(feasible,integer)
        print(row,flush=True);(out/'response_discovery.json').write_text(json.dumps(records,indent=2)+'\n')
    upper_results=[];best_upper=None
    for degree in upper_degrees:
        blocks=[];energy=[]
        for i,p in enumerate(P[:degree+1]):
            grows=[];hrows=[]
            for j,q in enumerate(P[:degree+1]):
                scale=12**(i+j)
                grows.append(np.array([[sum(a*b*M[k+l][r][s] for k,a in enumerate(p) for l,b in enumerate(q))/scale for s in range(14)] for r in range(14)]))
                hrows.append(np.array([[sum(a*b*M[k+l+1][r][s] for k,a in enumerate(p) for l,b in enumerate(q))/scale for s in range(14)] for r in range(14)]))
            blocks.append(grows);energy.append(hrows)
        G=np.block(blocks);H=np.block(energy);values,T=eigh(G);keep=values>max(values)*1e-11;T=T[:,keep]/np.sqrt(values[keep])
        ev,Y=eigh(T.T@H@T);weights=(T@Y[:,0]).reshape(degree+1,14)
        if chebyshev_upper:
            amplitude_scale=10**14/np.max(np.abs(weights))
            coefficients=[[round(float(x)*amplitude_scale) for x in row] for row in weights]
            upper,norm=polynomial_upper(M,upper_monomials(coefficients))
        else:
            coefficients=[[round(sum(F(float(weights[i,j]))*P[i][k]/12**i for i in range(k,degree+1))*10**14) for j in range(14)] for k in range(degree+1)]
            upper,norm=polynomial_upper(M,coefficients)
        row={'degree':degree,'upper':str(upper),'upper_float':float(upper),'numerical_rank':int(sum(keep)),'coefficient_count':14*(degree+1)}
        upper_results.append(row);print(row,flush=True)
        if best_upper is None or upper<best_upper[0]:best_upper=(upper,coefficients)
        (out/'upper_discovery.json').write_text(json.dumps(upper_results,indent=2)+'\n')
    if best is None:raise ValueError('No feasible response proposal')
    c=json.loads((ROOT/'singlet_energy_integer/certificate.json').read_text());c.pop('independent_upper')
    c.update(kind='h8_singlet_moment_energy_research_v1',lower=str(best[0]),response_polynomial=best[1])
    c['upper_chebyshev_coefficients' if chebyshev_upper else 'upper_polynomial_coefficients']=best_upper[1]
    (out/'candidate.json').write_text(json.dumps(c,indent=2)+'\n');receipt=replay(c)
    (out/'certificate.json').write_text(json.dumps(c,indent=2)+'\n');(out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print('EXACT ACCEPTED',receipt['lower'],receipt['upper_float'],'seconds',time.monotonic()-started,flush=True)
if __name__=='__main__':main()
