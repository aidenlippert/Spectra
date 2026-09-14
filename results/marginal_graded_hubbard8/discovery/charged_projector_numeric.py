"""Bounded discovery for a charged-sector extendibility constraint."""
from pathlib import Path
from fractions import Fraction as F
import sys, json, math, time
import numpy as np
from scipy.linalg import eigh
from scipy.optimize import minimize
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT))
from experiments.marginal_local_hubbard_block import sector_matrices
from experiments.marginal_charged_projectors import charged_vectors, charged_overlap_grams
BASE=ROOT/'results/marginal_graded_hubbard8';OUT=BASE/'charged_projector'


def main():
    start=time.monotonic();OUT.mkdir(exist_ok=True)
    base=json.loads((BASE/'six_site_projector/refined_certificate.json').read_text())
    local=base['local_window']
    sectors=sector_matrices(6,F(local['U']),F(local['t']),list(map(F,local['onsite_profile'])),
                            list(map(F,local['hopping_profile'])),F(local['V']),list(map(F,local['density_profile'])))
    _,K,cols=sectors[(2,3,1)]
    scale=np.sqrt([sum(a*a for a in col.values()) for col in cols])
    A=np.array(K,float)/scale[:,None]/scale[None,:]
    values,eigenvectors=eigh(A,subset_by_index=[0,0]);coeff=eigenvectors[:,0]/scale
    coeff=coeff/max(abs(coeff))
    vector={}
    for amplitude,col in zip(coeff,cols):
        z=round(float(amplitude)*10**8)
        for state,a in col.items():
            if z:vector[str(state)]=z*a
    (OUT/'source.json').write_text(json.dumps({'physical_state':vector,'source_sector':[2,3,1]},indent=2)+'\n')
    vectors,norm=charged_vectors(vector)
    ceilings=[]
    for m in (2,3,4):
        n,grams=charged_overlap_grams(vector,m)
        maximum=max(float(eigh(np.array(G,float)/float(n),subset_by_index=[len(G)-1,len(G)-1],eigvals_only=True)[0]) for G in grams.values())
        B=F(math.ceil((maximum+1e-8)*10**9),10**9)
        ceilings.append({'windows':m,'ceiling':str(B),'theta':str(B/m),'numerical_maximum':maximum,
                         'columns':sum(map(len,grams.values())),'maximum_block':max(map(len,grams.values()))})
    best=min(ceilings,key=lambda c:F(c['theta']));theta_ch=float(F(best['theta']))
    half={int(s):a for s,a in base['vector'].items()};half_norm=sum(a*a for a in half.values())
    matrices=[]
    for key,(_,K,cols) in sectors.items():
        scale=np.sqrt([sum(a*a for a in col.values()) for col in cols])
        A=np.array(K,float)/scale[:,None]/scale[None,:]
        def projector(v,n):
            w=np.array([sum(a*v.get(s,0) for s,a in col.items()) for col in cols],float)/scale/np.sqrt(float(n))
            return np.outer(w,w)
        matrices.append((key,A,projector(half,half_norm),sum((projector(v,norm) for v in vectors),np.zeros_like(A))))
    theta_half=float(F(base['projector_sum_ceiling'])/base['windows']);calls=0
    def minimum(k,l):
        return min((float(eigh(A+k*P+l*Q,subset_by_index=[0,0],eigvals_only=True)[0]),key) for key,A,P,Q in matrices)
    def objective(x):
        nonlocal calls
        calls+=1
        if calls>100:raise ValueError('Hard objective cap exceeded')
        k,l=x
        if min(x)<0 or max(x)>2:return 1000+sum(abs(x))
        return -(minimum(k,l)[0]-k*theta_half-l*theta_ch)/5
    x0=[float(F(base['penalty'])),0.02]
    result=minimize(objective,x0,method='Nelder-Mead',options={'maxfev':100,'xatol':1e-8,'fatol':1e-9})
    k,l=(F(round(float(x)*10**6),10**6) for x in result.x)
    ell,sector=minimum(float(k),float(l))
    lower=F(math.floor(ell*10**7)-1,10**7)
    candidate=dict(base,kind='hubbard_projector_extension_v3',penalty=str(k),penalized_lower=str(lower),
                   charged={'vector':vector,'windows':best['windows'],'projector_sum_ceiling':best['ceiling'],'penalty':str(l)})
    periodic=(lower-k*F(base['projector_sum_ceiling'])/base['windows']-l*F(best['theta']))/5
    old=(F(base['penalized_lower'])-F(base['penalty'])*F(base['projector_sum_ceiling'])/base['windows'])/5
    receipt={'ceilings':ceilings,'candidate_penalties':[str(k),str(l)],'numerical_local_minimum':ell,
             'active_sector':list(sector),'proposed_periodic_density':str(periodic),'previous_periodic_density':str(old),
             'improvement':str(periodic-old),'objective_evaluations':calls,'seconds':time.monotonic()-start,
             'scope':'Numerical fixed-profile discovery only; no all-state or ground-energy acceptance until exact replay.'}
    (OUT/'certificate.json').write_text(json.dumps(candidate,indent=2)+'\n')
    (OUT/'numeric_proposal.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt),flush=True)


if __name__=='__main__':main()
